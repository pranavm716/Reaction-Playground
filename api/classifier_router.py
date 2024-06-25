from fastapi import APIRouter, HTTPException

from api.utils import get_mol_and_image_encoding

from backend.computations import get_substructure_classifications

router = APIRouter(prefix="/classifier", tags=["classifier"])


@router.get("/")
def get_mol_classifications(smiles: str) -> list[str]:
    try:
        mol, _ = get_mol_and_image_encoding(smiles)
    except (ValueError, TypeError) as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    return get_substructure_classifications(mol)
