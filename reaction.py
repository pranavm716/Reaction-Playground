from dataclasses import dataclass

import rdkit.Chem.rdChemReactions as rd
from rdkit.Chem import AllChem


@dataclass(frozen=True, slots=True)
class Reaction:
    """
    A class that holds all the necessary information for a reaction.
    This class defines reactions more broadly than RDKit.
    While RDKit reactions precisely define which atoms and bonds are involved,
    this class organizes groups of specific reactions by behavior
    (ex. a more general 'oxidation' reaction that applies to many classes of
    molecules vs. any specific oxidation reaction). Here, the specific cases of the
    broader reaction are called subreactions.
    """

    name: str
    subreactions: list[rd.ChemicalReaction]
    description: str
    num_reactants: int

    @property
    def smarts_list(self):
        return [AllChem.ReactionToSmarts(rxn) for rxn in self.subreactions]

    def __str__(self):
        return f"{self.name}: {self.smarts_list}"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Reaction):
            return False

        return (
            self.name == other.name
            and set(self.smarts_list) == set(other.smarts_list)
            and self.description == other.description
            and self.num_reactants == other.num_reactants
        )
