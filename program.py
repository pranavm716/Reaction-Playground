from dataclasses import dataclass

from rdkit.Chem.rdchem import Mol

from reaction import Reaction
from ui import UI


@dataclass
class Program:
    ui: UI
    all_reactions: list[Reaction]
    start_mol: Mol
    target_mol: Mol | None = None

    def run(self):
        solver_mode = bool(self.target_mol)
