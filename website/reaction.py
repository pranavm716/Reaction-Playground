import pathlib
from functools import cached_property
from typing import Iterator

import rdkit.Chem.rdChemReactions as rd
import yaml
from pydantic import BaseModel, Field, model_validator, TypeAdapter, ConfigDict
from rdkit.Chem import AllChem


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
    smarts_list: list[str]
    subreactions: list[rd.ChemicalReaction] = Field(
        description="A list of the more specific rdkit reactions that are sub-scenarios "
        "of this more general reaction object."
    )
    multiple_reactants_prompts: list[str] | None = Field(
        default=None,
        description="<...> Should only be set for reactions with more than one reactant.",  # TODO
    )
    description: str

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


def read_all_reactions_from_file(path: pathlib.Path) -> list[Reaction]:
    """
    Returns a list of all available reactions.
    """

    reaction_list = yaml.safe_load(path.read_text())
    ta = TypeAdapter(list[Reaction])
    all_reactions = ta.validate_python(reaction_list)
    return all_reactions
