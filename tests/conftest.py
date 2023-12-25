import pytest

from website.config import ALL_REACTIONS_FILE_PATH, DISABLE_RDKIT_WARNINGS
from website.reaction import Reaction, read_all_reactions_from_file
import pathlib


@pytest.fixture
def all_reactions() -> list[Reaction]:
    all_reactions_file_path = (
        pathlib.Path(__file__).resolve().parent.parent
        / "website"
        / ALL_REACTIONS_FILE_PATH
    )
    real_reactions = read_all_reactions_from_file(all_reactions_file_path)

    made_up_reaction = {
        "name": "Made up reaction",
        "subreactions": [
            {
                "reactants": [
                    {"smarts": "[CX3H1](=[O])[#6]", "classification": "something"},
                    {"smarts": "[#6][CH]=[CH2]", "classification": "something"},
                    {"smarts": "[CH,CH2][OH]", "classification": "something"},
                ],
                "products": [{"smarts": "[C]", "classification": "something"}],
            }
        ],
        "description": "Some description",
    }

    return real_reactions + [Reaction(**made_up_reaction)]


@pytest.fixture(autouse=True)
def disable_rdkit_warnings_if_configured() -> None:
    if DISABLE_RDKIT_WARNINGS:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.warning")
