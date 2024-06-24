from fastapi import APIRouter

from api.utils import get_mol_and_image_encoding

from backend.computations import get_substructure_classifications

router = APIRouter(prefix="/classifier", tags=["classifier"])


@router.get("/")
def get_mol_classifications(mol_smiles: str) -> list[str]:
    mol, _ = get_mol_and_image_encoding(mol_smiles)
    return get_substructure_classifications(mol)
