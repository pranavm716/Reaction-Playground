import pathlib

import yaml

import rdkit.Chem.rdChemReactions as rd
from rdkit.Chem import AllChem
from pydantic import BaseModel, Field, model_validator, TypeAdapter, ConfigDict
from functools import cached_property

from website.config import DISABLE_RDKIT_WARNINGS

BASE_REACTION_MODEL_CONFIG = ConfigDict(frozen=True, arbitrary_types_allowed=True)


class Subreaction(BaseModel):
    model_config = BASE_REACTION_MODEL_CONFIG

    reaction: rd.ChemicalReaction = Field(
        description="The reaction object representing this subreaction, created by parsing the SMARTS string."
    )
    smarts: str = Field(description="The SMARTS string representing this reaction.")
    effect: str = Field(
        description="The effect of this specific rdkit defined reaction.",
        examples=["Hydrolyze the amide into a carboxylic acid and an amine."],
    )

    @model_validator(mode="before")
    def transform_smarts_string(cls, values):
        smarts_reaction = AllChem.ReactionFromSmarts(values["smarts"])
        smarts_reaction.Initialize()
        values["reaction"] = smarts_reaction
        return values


class Reaction(BaseModel):
    """
    A class that holds all the necessary information for a reaction.
    This class defines reactions more broadly than RDKit.
    While RDKit reactions precisely define which atoms and bonds are involved,
    this class organizes groups of specific reactions by behavior
    (ex. a more general 'oxidation' reaction that applies to many classes of
    molecules vs. any specific oxidation reaction). Here, the specific cases of the
    broader reaction are called subreactions.
    """

    model_config = BASE_REACTION_MODEL_CONFIG

    name: str
    subreactions: list[Subreaction] = Field(
        description="A list of the more specific rdkit reactions that are sub-scenarios of this more general reaction object."
    )
    multiple_reactants_prompts: list[str] | None = Field(
        default=None,
        description="<...> Should only be set for reactions with more than one reactant.",  # TODO
    )
    description: str

    @cached_property
    def smarts_list(self) -> list[str]:
        return [subreaction.smarts for subreaction in self.subreactions]

    @cached_property
    def num_reactants(self) -> int:
        # All the subreactions have the same number of reactants, so we can just use the first one
        return self.subreactions[0].reaction.GetNumReactantTemplates()  # noqa

    def __str__(self):
        return f"{self.name}: {self.smarts_list}"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Reaction) and set(self.smarts_list) == set(
            other.smarts_list
        )


def read_all_reactions_from_file(path: pathlib.Path) -> list[Reaction]:
    """
    Returns a list of all available reactions.
    """
    if DISABLE_RDKIT_WARNINGS:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.warning")

    reaction_dict = yaml.safe_load(path.read_text())
    ta = TypeAdapter(list[Reaction])
    all_reactions = ta.validate_python(reaction_dict)
    return all_reactions
