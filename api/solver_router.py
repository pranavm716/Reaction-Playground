from fastapi import APIRouter, HTTPException, Query

from api.utils import get_mol_and_image_encoding, mol_to_base64
from pydantic import BaseModel

from backend.computations import (
    ALL_REACTIONS,
    copy_mol,
    find_synthetic_pathway,
    generate_multi_step_product,
)
from rdkit import Chem


router = APIRouter(prefix="/solver", tags=["solver"])


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
        raise HTTPException(status_code=400, detail=str(exc))

    return encoding


# Mapping of the step number to a list of image encodings of the products generated
# by running the corresponding reaction for that step
SolverModeImageData = dict[int, list[str]]


class SolverModeResponse(BaseModel):
    path_found: bool
    num_steps: int
    reaction_names: list[str]
    choice_pathway: list[int]
    solver_image_encodings: SolverModeImageData


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
    reaction_names: list[str] = []
    current_mol = copy_mol(start_mol)
    for step_number, (reaction_key, choice) in enumerate(
        zip(reaction_pathway, choice_pathway)
    ):
        reaction_names.append(ALL_REACTIONS[reaction_key].name)
        products = generate_multi_step_product(current_mol, reaction_key)

        product_images = [mol_to_base64(product) for product in products]
        solver_image_encodings[step_number] = product_images

        current_mol = products[choice]

    return SolverModeResponse(
        path_found=path_found,
        num_steps=num_steps,
        reaction_names=reaction_names,
        choice_pathway=choice_pathway,
        solver_image_encodings=solver_image_encodings,
    )
