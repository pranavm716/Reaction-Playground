from fastapi import APIRouter, HTTPException, Query

from api.utils import get_mol_and_image_encoding


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
