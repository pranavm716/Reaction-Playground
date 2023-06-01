from abc import ABC, abstractmethod

from config import Mol2dTuple
from rdkit.Chem.rdchem import Mol



class UI(ABC):
    """A generic UI abstract class."""

    @abstractmethod
    def prompt_for_start_and_target_molecule(self) -> tuple[Mol, Mol]:
        """
        Prompts the user for the SMILES of the starting and target molecule.
        The target molecule is optional.
        """
        pass

    @abstractmethod
    def display_mol(self, mol: Mol) -> None:
        """
        Draws a single molecule on the screen.
        """
        pass

    @abstractmethod
    def display_2d_mol_tuple(self, products: Mol2dTuple) -> None:
        """
        Draws a grid of molecules on screen.
        """
        pass
