from fastapi import FastAPI
from rdkit import Chem
from urllib.parse import quote

from rdkit.Chem.rdchem import Mol


def start_and_target_mols_are_valid(
    start_mol_smiles: str, target_mol_smiles: str
) -> bool:
    try:
        start_mol = Chem.MolFromSmiles(start_mol_smiles)
        Chem.SanitizeMol(start_mol)

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


def get_mol_from_smiles(smiles: str) -> Mol:
    return Chem.MolFromSmiles(smiles)
