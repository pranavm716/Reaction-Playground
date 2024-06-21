from rdkit.Chem import Draw
from rdkit.Chem.rdchem import Mol
from io import BytesIO
import base64
from PIL.Image import Image as PILImage


def mol_to_base64(mol: Mol) -> str:
    return _img_to_base64(_construct_mol_image(mol))


def _img_to_base64(img: PILImage) -> str:
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    image_data = base64.b64encode(buffered.getvalue()).decode("utf-8")
    return image_data


def _construct_mol_image(mol: Mol) -> PILImage:
    return Draw.MolToImage(mol)


# old
# def start_and_target_mols_are_valid(
#     start_mol_smiles: str, target_mol_smiles: str | None = None
# ) -> bool:
#     try:
#         start_mol = Chem.MolFromSmiles(start_mol_smiles)
#         Chem.SanitizeMol(start_mol)

#         if target_mol_smiles:
#             target_mol = Chem.MolFromSmiles(target_mol_smiles)
#             Chem.SanitizeMol(target_mol)
#     except (TypeError, ValueError):
#         return False
#     else:
#         return True
