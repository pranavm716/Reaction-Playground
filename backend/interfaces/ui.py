from abc import ABC, abstractmethod

from rdkit.Chem.rdchem import Mol

from backend.datatypes import Mol2dTuple, MolTuple
from backend.reaction import Reaction


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
    def display_mol_tuple(self, mols: MolTuple) -> None:
        """
        Displays a 1d tuple of molecules in a line.
        """
        pass

    @abstractmethod
    def get_user_input(self, prompt: str = "") -> str:
        pass

    @abstractmethod
    def print(self, text: str = "") -> None:
        pass

    # -------------------- For solver mode --------------------
    @abstractmethod
    def display_solver_mode_intro(self, start_mol: Mol, target_mol: Mol) -> None:
        """
        Displays the solver mode introductory text and images.
        """
        pass

    @abstractmethod
    def display_solver_step(
        self,
        current_mol: Mol,
        step_number: int,
        reaction: Reaction,
        products: MolTuple,
        choice: int,
    ) -> None:
        """
        Displays information about the current step of the solver.
        Displays which reaction turns the current mol into the products and
        which product to pick for the next step.
        """
        pass

    # ---------------------------------------------------------

    # -------------------- For playground mode ----------------
    @abstractmethod
    def display_2d_mol_tuple(self, products: Mol2dTuple) -> None:
        """
        Draws a 2d grid of molecules on screen.
        """
        pass

    @abstractmethod
    def display_compatible_reactions(
        self, start_mol: Mol, possible_reactions: list[Reaction]
    ) -> None:
        """
        Displays the start molecule and the list of compatible reactions.
        """
        pass

    @abstractmethod
    def prompt_for_scenario_index(self, products: Mol2dTuple) -> int:
        """
        When multiple scenarios are generated from the results of single step mode,
        this method prompts the user for the scenario index.
        """
        pass

    @abstractmethod
    def prompt_for_product_index(self, scenario: MolTuple) -> int:
        """
        For a given scenario, this method prompts the user for the product index.
        """
        pass

    # ---------------------------------------------------------
