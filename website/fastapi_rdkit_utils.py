from fastapi import FastAPI
from rdkit import Chem
from urllib.parse import quote

from rdkit.Chem import Draw
from rdkit.Chem.rdchem import Mol
from io import BytesIO
import base64
from PIL.Image import Image as PILImage


def start_and_target_mols_are_valid(
    start_mol_smiles: str, target_mol_smiles: str | None = None
) -> bool:
    try:
        start_mol = Chem.MolFromSmiles(start_mol_smiles)
        Chem.SanitizeMol(start_mol)

        if target_mol_smiles:
            target_mol = Chem.MolFromSmiles(target_mol_smiles)
            Chem.SanitizeMol(target_mol)
    except (TypeError, ValueError):
        return False
    else:
        return True


def construct_query_url(app: FastAPI, url_path_for: str, **query_params: str) -> str:
    url = app.url_path_for(url_path_for)
    query_string = [f"{k}={quote(v)}" for k, v in query_params.items()]
    return f"{url}?{'&'.join(query_string)}"


def img_to_base64(img: PILImage) -> str:
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    image_data = base64.b64encode(buffered.getvalue()).decode("utf-8")
    return image_data


def construct_mol_image(mol: Mol) -> PILImage:
    return Draw.MolToImage(mol)


def mol_to_base64(mol: Mol) -> str:
    return img_to_base64(construct_mol_image(mol))
