import os

import pytest
from rdkit.Chem import AllChem

from program import read_all_reactions_from_file
from reaction import Reaction


@pytest.fixture
def all_reactions() -> list[Reaction]:
    all_reactions_file_path = os.path.join(
        os.path.dirname(__file__), "../all_reactions.txt"
    )
    real_reactions = read_all_reactions_from_file(all_reactions_file_path)
    made_up_reaction = Reaction(
        name="Made up reaction",
        subreactions=[
            AllChem.ReactionFromSmarts(
                "[CX3H1:1](=[O:2])[#6:3].[#6:1][CH:2]=[CH2:3].[CH,CH2:1][OH:2]>>[Au]"
            )
        ],
        description="Some description",
        num_reactants=3,
    )
    return real_reactions + [made_up_reaction]
