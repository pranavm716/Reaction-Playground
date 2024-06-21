from fastapi import APIRouter, HTTPException

from backend.computations import find_possible_reaction_keys
from backend.fastapi_rdkit_utils import smiles_to_base64
from backend.reaction import ReactionKey
from rdkit import Chem


router = APIRouter(prefix="/playground", tags=["playground"])


@router.get("/step")
def get_mol_image_and_valid_reactions(smiles: str) -> tuple[str, list[ReactionKey]]:
    try:
        encoding = smiles_to_base64(smiles)
    except (
        TypeError,
        ValueError,
    ) as exc:  # Something went wrong while parsing the SMILES
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    reaction_keys = find_possible_reaction_keys(
        Chem.MolFromSmiles(smiles), solver_mode=False
    )

    return encoding, reaction_keys
