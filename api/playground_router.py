from fastapi import APIRouter, HTTPException

from backend.computations import find_possible_reaction_keys, get_reactions_from_keys
from backend.fastapi_rdkit_utils import smiles_to_base64
from backend.reaction import Reaction
from rdkit import Chem


router = APIRouter(prefix="/playground", tags=["playground"])


@router.get("/step")
def get_mol_image_and_valid_reactions(smiles: str) -> tuple[str, list[Reaction]]:
    try:
        encoding = smiles_to_base64(smiles)
    except (TypeError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    reaction_keys = find_possible_reaction_keys(
        Chem.MolFromSmiles(smiles), solver_mode=False
    )
    reactions = get_reactions_from_keys(reaction_keys)

    return encoding, reactions
