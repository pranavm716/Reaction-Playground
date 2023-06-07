import os

import pytest
from rdkit.Chem import AllChem

from backend.config import Config, read_config
from backend.program import read_all_reactions_from_file
from backend.reaction import Reaction


@pytest.fixture
def all_reactions() -> list[Reaction]:
    all_reactions_file_path = os.path.join(
        os.path.dirname(__file__), "../all_reactions.txt"
    )
    real_reactions = read_all_reactions_from_file(all_reactions_file_path)

    made_up_subreaction = AllChem.ReactionFromSmarts(
        "[CX3H1](=[O])[#6].[#6][CH]=[CH2].[CH,CH2][OH]>>[C]"
    )
    made_up_subreaction.Initialize()
    made_up_reaction = Reaction(
        name="Made up reaction",
        subreactions=[made_up_subreaction],
        description="Some description",
        num_reactants=3,
    )
    return real_reactions + [made_up_reaction]


@pytest.fixture
def config() -> Config:
    config_file_path = os.path.join(os.path.dirname(__file__), "../config.json")
    c = read_config(config_file_path)
    if c.disable_rdkit_warnings:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.warning")

    c.all_reactions_file_path = os.path.join(
        os.path.dirname(__file__), "..", c.all_reactions_file_path
    )
    return c
