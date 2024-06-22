import copy
import functools
from fastapi import APIRouter, HTTPException, Query, Body

from backend.computations import (
    ALL_REACTIONS,
    find_possible_reaction_keys,
    generate_multi_step_product,
    generate_single_step_product,
    get_reactant_position_of_mol_in_reaction,
    get_reactions_from_keys,
)
from backend.fastapi_rdkit_utils import mol_to_base64
from backend.reaction import Reaction, ReactionKey
from rdkit import Chem
from rdkit.Chem.rdchem import Mol


router = APIRouter(prefix="/playground", tags=["playground"])


@router.get("/start")
def get_mol_image_and_valid_reactions(
    smiles: str = Query(...),
) -> tuple[str, list[Reaction]]:
    """
    Takes a SMILES string and returns a 2-tuple containing:
    - the base64 encoding of the molecule's image
    - a list of reaction keys that can be performed on the molecule

    The SMILES string must contain only one valid molecule. An exception is raised otherwise.
    """

    # Validation
    if not smiles:
        raise HTTPException(status_code=400, detail="Drawing cannot be empty.")
    if "." in smiles:
        raise HTTPException(
            status_code=400,
            detail=f"Drawing must contain only one molecule. Received {smiles.count(".") + 1} molecules.",
        )

    mol = Chem.MolFromSmiles(smiles)
    encoding = mol_to_base64(mol)

    reaction_keys = find_possible_reaction_keys(mol, solver_mode=False)
    reactions = get_reactions_from_keys(reaction_keys)

    return encoding, reactions


@router.get("/missing-reactants")
def get_missing_reactant_prompts(
    smiles: str = Query(...),
    reaction_key: ReactionKey = Query(...),
) -> list[str]:
    """
    Returns the prompts for the given reaction that requires additional reactants.

    Pre: The reaction key is valid for the given molecule.
    Pre: The reaction key is for a reaction with multiple reactants.
    """
    reaction = ALL_REACTIONS[reaction_key]
    if reaction.multiple_reactants_prompts is None:
        raise HTTPException(
            status_code=400, detail="This reaction has only one reactant."
        )

    mol = Chem.MolFromSmiles(smiles)
    reactant_position = get_reactant_position_of_mol_in_reaction(mol, reaction_key)

    multiple_reactant_prompts = copy.copy(reaction.multiple_reactants_prompts)
    multiple_reactant_prompts.pop(reactant_position)

    return multiple_reactant_prompts


@router.get("/reaction/single-reactant")
def get_products_single_reactant(
    smiles: str = Query(...),
    reaction_key: ReactionKey = Query(...),
) -> tuple[tuple[str, str], ...]:
    """
    Takes a SMILES string and a reaction key and returns an arbitrary length products tuple containing 2-tuples of:
    - the base64 encoding of the product's image
    - the SMILES of the product
    Pre: The reaction key is valid for the given molecule.
    """
    products = generate_multi_step_product(Chem.MolFromSmiles(smiles), reaction_key)
    return tuple((mol_to_base64(p), Chem.MolToSmiles(p)) for p in products)


@router.post("/reaction/multiple-reactants")
def get_products_multiple_reactants(
    smiles: str = Body(...),
    extra_reactant_smiles: list[str] = Body(...),
    reaction_key: ReactionKey = Body(...),
):
    mol = Chem.MolFromSmiles(smiles)
    reactants: list[Mol] = [
        Chem.MolFromSmiles(smiles) for smiles in extra_reactant_smiles
    ]

    reactant_position = get_reactant_position_of_mol_in_reaction(mol, reaction_key)
    reactants.insert(reactant_position, mol)

    products = generate_single_step_product(tuple(reactants), reaction_key)
    if not products:
        raise HTTPException(
            status_code=400,
            detail="Invalid reactant molecules provided for this reaction.",
        )

    # NOTE: I (for now) merge a m x n 2d tuple of single step products into a 1 x m*n 1d tuple of products
    # I can't figure out how to implement/ best represent multi step product for multiple reactants
    products = functools.reduce(lambda x, y: x + y, products)

    return tuple((mol_to_base64(p), Chem.MolToSmiles(p)) for p in products)
