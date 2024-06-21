import pathlib
from enum import StrEnum, auto
from functools import cached_property
from typing import Iterator

import rdkit.Chem.rdChemReactions as rd
import yaml
from pydantic import BaseModel, Field, model_validator, ConfigDict
from rdkit.Chem import AllChem


class ReactionKey(StrEnum):
    """
    Enum class representing a reaction. Used internally only.
    """

    # Reactions involving a single reactant
    hydrolysis = auto()
    dmp_pcc_oxidation = auto()
    h2cro4_oxidation = auto()
    nabh4_reduction = auto()
    lialh4_reduction = auto()
    pbr3_bromination = auto()
    ozonolysis = auto()
    nacn_nitrile_synthesis = auto()
    om_dm = auto()
    hydroboration_oxidation = auto()
    grignard_reagent = auto()
    socl2_acid_chloride_synthesis = auto()
    deprotonation = auto()

    # Reactions involving multiple reactants
    grignard_reaction = auto()
    amide_synthesis_from_acid_chloride = auto()
    ester_synthesis_from_acid_chloride = auto()
    fischer_esterification = auto()
    williamson_ether_synthesis = auto()


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

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    name: str
    smarts_list: list[str] = Field(
        description="A list of the SMARTS strings that represent the subreactions "
        "for this reaction."
    )
    subreactions: list[rd.ChemicalReaction] = Field(
        exclude=True,
        description="A list of the more specific rdkit reactions that are sub-scenarios "
        "of this more general reaction object.",
    )
    multiple_reactants_prompts: list[str] | None = Field(
        default=None,
        description="List of prompts that ask the user to enter the SMILES of the missing reactants. "
        "Should only be set for reactions with more than one reactant.",
        examples=[
            "Enter the SMILES of carbon dioxide (O=C=O) or an aldehyde, ketone, or nitrile: "
        ],
    )
    description: str
    reaction_key: ReactionKey

    @model_validator(mode="before")
    def transform_smarts_list(cls, values):
        subreactions: list[rd.ChemicalReaction] = []
        for smarts in values["smarts_list"]:
            subreaction = AllChem.ReactionFromSmarts(smarts)
            subreaction.Initialize()
            subreactions.append(subreaction)
        values["subreactions"] = subreactions
        return values

    @cached_property
    def num_reactants(self) -> int:
        # All the subreactions have the same number of reactants, so we can just use the first one
        return self.subreactions[0].GetNumReactantTemplates()  # noqa

    def __iter__(self) -> Iterator[rd.ChemicalReaction]:
        """
        Makes it more convenient to iterate over the subreactions by yielding the underlying rd.ChemicalReaction objects
        """
        yield from self.subreactions


# Mapping of the reaction key to its pydantic Reaction object
ReactionDict = dict[ReactionKey, Reaction]


def read_all_reactions_from_file(path: pathlib.Path) -> ReactionDict:
    """
    Returns a dict of all available reactions.
    """

    data = yaml.safe_load(path.read_text())
    all_reactions: ReactionDict = {}
    for key, reaction_data in data.items():
        key = ReactionKey(key)
        all_reactions[key] = Reaction(**reaction_data, reaction_key=key)

    return all_reactions
