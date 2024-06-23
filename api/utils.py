from rdkit.Chem import Draw
from rdkit.Chem.rdchem import Mol
from io import BytesIO
import base64
from PIL.Image import Image as PILImage
from rdkit import Chem


def get_mol_and_image_encoding(smiles: str) -> tuple[Mol, str]:
    # Validation
    if not smiles:
        raise ValueError("Drawing cannot be empty.")
    if "." in smiles:
        raise ValueError(
            f"Drawing must contain only one molecule. Received {smiles.count(".") + 1} molecules.",
        )

    mol = Chem.MolFromSmiles(smiles)
    Chem.SanitizeMol(mol)
    encoding = mol_to_base64(mol)
    return mol, encoding


def mol_to_base64(mol: Mol) -> str:
    return _img_to_base64(_construct_mol_image(mol))


def _img_to_base64(img: PILImage) -> str:
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    image_data = base64.b64encode(buffered.getvalue()).decode("utf-8")
    return image_data


def _construct_mol_image(mol: Mol) -> PILImage:
    return Draw.MolToImage(mol)