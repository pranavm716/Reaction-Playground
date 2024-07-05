from fastapi import APIRouter, HTTPException, Query, Depends

from api.response import MolImageMetadata, SolverModeImageData, SolverModeResponse
from api.utils import get_mol_and_image_encoding, mol_to_base64

from backend.computations import (
    ALL_REACTIONS,
    copy_mol,
    find_synthetic_pathway,
    generate_multi_step_product,
)
from api.auth import verify_api_key
from backend.reaction import Reaction
from rdkit import Chem


router = APIRouter(
    prefix="/solver", tags=["solver"], dependencies=[Depends(verify_api_key)]
)


@router.get("/get-mol-image")
def get_image_encoding(
    smiles: str = Query(...),
) -> str:
    """
    Takes a SMILES string and returns the base64 encoding of the molecule's image.
    """
    try:
        _, encoding = get_mol_and_image_encoding(smiles)
    except (ValueError, TypeError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    return encoding


@router.get("/run")
def run_solver_mode(start_smiles: str, target_smiles: str) -> SolverModeResponse:
    """
    Runs the auto synthetic pathway solver. Returns information about the calculated synthetic pathway,
    including encodings of images of the molecules at each step.

    Pre: start_smiles and target_smiles represent valid molecules.
    """

    start_mol = Chem.MolFromSmiles(start_smiles)
    target_mol = Chem.MolFromSmiles(target_smiles)

    # Includes the target molecule as the last entry, does not include the starting molecule
    solver_image_encodings: SolverModeImageData = {}

    path_found, reaction_pathway, choice_pathway = find_synthetic_pathway(
        start_mol, target_mol
    )
    num_steps = len(choice_pathway)

    # Reconstruct the synthetic pathway
    reactions: list[Reaction] = []
    current_mol = copy_mol(start_mol)
    for step_number, (reaction_key, choice) in enumerate(
        zip(reaction_pathway, choice_pathway)
    ):
        reactions.append(ALL_REACTIONS[reaction_key])
        products = generate_multi_step_product(current_mol, reaction_key)

        product_images: list[MolImageMetadata] = [
            MolImageMetadata(
                smiles=Chem.MolToSmiles(product),
                encoding=mol_to_base64(product),
            )
            for product in products
        ]
        solver_image_encodings[step_number] = product_images

        current_mol = products[choice]

    return SolverModeResponse(
        path_found=path_found,
        num_steps=num_steps,
        reactions=reactions,
        choice_pathway=choice_pathway,
        starting_encoding=mol_to_base64(start_mol),
        target_encoding=mol_to_base64(target_mol),
        solver_image_metadata=solver_image_encodings,
    )
