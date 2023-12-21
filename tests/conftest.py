import pytest

from website.config import ALL_REACTIONS_FILE_PATH
from website.reaction import Reaction, read_all_reactions_from_file, Subreaction
import pathlib


@pytest.fixture
def all_reactions() -> list[Reaction]:
    all_reactions_file_path = (
        pathlib.Path(__file__).resolve().parent.parent
        / "website"
        / ALL_REACTIONS_FILE_PATH
    )
    real_reactions = read_all_reactions_from_file(all_reactions_file_path)

    made_up_reaction = Reaction(
        name="Made up reaction",
        subreactions=[
            Subreaction(
                smarts="[CX3H1](=[O])[#6].[#6][CH]=[CH2].[CH,CH2][OH]>>[C]",
                effect="A made up reaction.",
            )
        ],
        description="Some description",
    )
    return real_reactions + [made_up_reaction]
