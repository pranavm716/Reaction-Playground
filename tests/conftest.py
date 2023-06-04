import os

import pytest

from program import read_all_reactions_from_file
from reaction import Reaction


@pytest.fixture
def all_reactions() -> list[Reaction]:
    all_reactions_file_path = os.path.join(
        os.path.dirname(__file__), "../all_reactions.txt"
    )
    return read_all_reactions_from_file(all_reactions_file_path)
