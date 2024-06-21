from fastapi import APIRouter, HTTPException

from backend.computations import (
    find_possible_reaction_keys,
    generate_multi_step_product,
    get_reactions_from_keys,
)
from backend.fastapi_rdkit_utils import mol_to_base64
from backend.reaction import Reaction, ReactionKey
from rdkit import Chem


router = APIRouter(prefix="/playground", tags=["playground"])


@router.get("/step-start")
def get_mol_image_and_valid_reactions(smiles: str) -> tuple[str, list[Reaction]]:
    """
    Takes a SMILES string and returns a 2-tuple containing:
    - the base64 encoding of the molecule's image
    - a list of reaction keys that can be performed on the molecule

    The SMILES string must contain only one valid molecule. An exception is raised otherwise.
    """

    # Validation
    if not smiles:
        raise HTTPException(status_code=400, detail="SMILES cannot be empty.")
    if "." in smiles:
        raise HTTPException(
            status_code=400,
            detail=f"SMILES must contain only one molecule. Received {smiles.count(".") + 1} molecules.",
        )

    mol = Chem.MolFromSmiles(smiles)
    encoding = mol_to_base64(mol)

    reaction_keys = find_possible_reaction_keys(mol, solver_mode=False)
    reactions = get_reactions_from_keys(reaction_keys)

    return encoding, reactions


@router.get("/step-reaction")
def get_products(smiles: str, reaction_key: ReactionKey) -> tuple[tuple[str, str], ...]:
    """
    Takes a SMILES string and a reaction key and returns an arbitrary length products tuple containing 2-tuples of:
    - the base64 encoding of the product's image
    - the SMILES of the product
    Pre: The reaction key is valid for the given molecule.
    """
    products = generate_multi_step_product(Chem.MolFromSmiles(smiles), reaction_key)
    return tuple((mol_to_base64(p), Chem.MolToSmiles(p)) for p in products)
