import pathlib
from functools import cached_property
from typing import Iterator

import rdkit.Chem.rdChemReactions as rd
import yaml
from pydantic import BaseModel, Field, model_validator, TypeAdapter, ConfigDict
from rdkit.Chem import AllChem


class Substructure(BaseModel):
    model_config = ConfigDict(frozen=True)

    smarts: str
    classification: str


class Subreaction(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    reactants: list[Substructure]
    products: list[Substructure]
    reaction: rd.ChemicalReaction = Field(
        description="The reaction object representing this subreaction, created by parsing the SMARTS string."
    )

    @model_validator(mode="before")
    def transform_smarts_string(cls, values):
        reactant_smarts = ".".join(r["smarts"] for r in values["reactants"])
        product_smarts = ".".join(p["smarts"] for p in values["products"])
        subreaction = AllChem.ReactionFromSmarts(
            ">>".join((reactant_smarts, product_smarts))
        )
        subreaction.Initialize()
        values["reaction"] = subreaction
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

    model_config = ConfigDict(frozen=True)

    name: str
    subreactions: list[Subreaction] = Field(
        description="A list of the more specific rdkit reactions that are sub-scenarios "
        "of this more general reaction object."
    )
    multiple_reactants_prompts: list[str] | None = Field(
        default=None,
        description="<...> Should only be set for reactions with more than one reactant.",  # TODO
    )
    description: str

    @cached_property
    def smarts_list(self) -> list[str]:
        return [
            AllChem.ReactionToSmarts(subreaction.reaction)
            for subreaction in self.subreactions
        ]

    @cached_property
    def num_reactants(self) -> int:
        # All the subreactions have the same number of reactants, so we can just use the first one
        return len(self.subreactions[0].reactants)

    def __iter__(self) -> Iterator[rd.ChemicalReaction]:
        """
        Makes it more convenient to iterate over the subreactions by yielding the underlying rd.ChemicalReaction objects
        """
        yield from (subreaction.reaction for subreaction in self.subreactions)

    def __str__(self):
        return f"{self.name}: {self.smarts_list}"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Reaction) and set(self.smarts_list) == set(
            other.smarts_list
        )

    def __hash__(self):
        return hash(tuple(self.smarts_list))


def read_all_reactions_from_file(path: pathlib.Path) -> list[Reaction]:
    """
    Returns a list of all available reactions.
    """

    reaction_list = yaml.safe_load(path.read_text())
    ta = TypeAdapter(list[Reaction])
    all_reactions = ta.validate_python(reaction_list)
    return all_reactions
