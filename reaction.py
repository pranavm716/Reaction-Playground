from dataclasses import dataclass

from rdkit.Chem import AllChem
import rdkit.Chem.rdChemReactions as rd


@dataclass(frozen=True, slots=True)
class Reaction:
    """
    A class that holds all the necessary information for a reaction.
    This class defines reactions more broadly than RDKit.
    While RDKit reactions precisely define which atoms and bonds are involved,
    this class organizes groups of specific reactions by behavior
    (ex. a more general 'oxidation' reaction that applies to many classes of 
    molecules vs. any specific oxidation reaction).
    """

    name: str
    reactions_list: list[rd.ChemicalReaction]
    description: str
    num_reactants: int

    @property
    def smarts_list(self):
        return [AllChem.ReactionToSmarts(rxn) for rxn in self.reactions_list]

    def __str__(self):
        return f"{self.name}: {self.smarts_list}"
