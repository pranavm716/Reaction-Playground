from abc import ABC, abstractmethod

from IPython.display import display
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdchem import Mol

from datatypes import Mol2dTuple, MolTuple


class UI(ABC):
    """An abstract class representing the user interface for this project."""

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
    def display_solver_mode_intro(self, start_mol: Mol, target_mol: Mol) -> None:
        """
        Displays the solver mode introductory text and images.
        """
        pass

    @abstractmethod
    def display_2d_mol_tuple(self, products: Mol2dTuple) -> None:
        """
        Draws a 2d grid of molecules on screen.
        """
        pass


class GoogleColabUI(UI):
    """A concrete implementation of the UI class. This class contains methods
    to gather user input and display text/ images in the Google Colab terminal.
    """

    def prompt_for_start_and_target_molecule(self) -> tuple[Mol, Mol | None]:
        start_smiles = input("Enter the SMILES of the starting molecule: ")
        start_mol = Chem.MolFromSmiles(start_smiles)

        target_mol = None
        if input("Would you like to input a target molecule? (y/n): ") == "y":
            target_smiles = input("Enter the SMILES of the target molecule: ")
            target_mol = Chem.MolFromSmiles(target_smiles)

        return start_mol, target_mol

    def display_mol(self, mol: Mol) -> None:
        display(Draw.MolToImage(mol))

    def display_solver_mode_intro(self, start_mol: Mol, target_mol: Mol) -> None:
        print(
            "\nThe goal is to find a reaction pathway that converts the molecule on the left into the molecule on the right."
        )
        display(Draw.MolsToGridImage([start_mol, target_mol]))

    def display_2d_mol_tuple(self, products: Mol2dTuple) -> None:
        for index, scenario in enumerate(products, start=1):
            print(f"Scenario #{index}")
            display(Draw.MolsToGridImage(list(scenario), molsPerRow=len(scenario)))
            print("\n")
