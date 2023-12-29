import pathlib

import pytest

from website.config import ALL_REACTIONS_FILE_PATH, ALL_SUBSTRUCTURES_FILE_PATH
from website.datatypes import SubstructureDict
from website.reaction import (
    Reaction,
    read_all_reactions_from_file,
    read_all_substructures_from_file,
)


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
        "smarts_list": [
            "[CX3H1:1](=[O])[#6].[#6:2][CH]=[CH2].[CH,CH2:3][OH]>>[C:1].[#6:2].[C:3]"
        ],
        "description": "Some description",
    }

    return real_reactions + [Reaction(**made_up_reaction)]


@pytest.fixture
def all_substructures() -> SubstructureDict:
    substructures_file_path = (
        pathlib.Path(__file__).resolve().parent.parent
        / "website"
        / ALL_SUBSTRUCTURES_FILE_PATH
    )
    return read_all_substructures_from_file(substructures_file_path)
