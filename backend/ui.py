from abc import ABC, abstractmethod

from rdkit import Chem
from rdkit.Chem import Draw
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


class CLI(UI):
    """A concrete implementation of the UI class. This class contains methods
    to gather user input, display images, and display text in the terminal.
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
        img = Draw.MolToImage(mol)
        # display(img)
        img.show()

    def display_mol_tuple(self, mols: MolTuple) -> None:
        img = Draw.MolsToGridImage(list(mols), molsPerRow=len(mols))
        # display(img)
        img.show()

    def get_user_input(self, prompt: str = "") -> str:
        return input(prompt)

    def print(self, text: str = "") -> None:
        print(text)

    # -------------------- For solver mode ----------------
    def display_solver_mode_intro(self, start_mol: Mol, target_mol: Mol) -> None:
        print(
            "\nThe goal is to find a reaction pathway that converts the molecule on the left into the molecule on the right."
        )
        self.display_mol_tuple((start_mol, target_mol))

    def display_solver_step(
        self,
        current_mol: Mol,
        step_number: int,
        reaction: Reaction,
        products: MolTuple,
        choice: int,
    ) -> None:
        self.display_mol(current_mol)
        print(f"\nStep #{step_number} - {reaction.name}.", end=" ")

        num_products = len(products)
        if num_products == 1:
            print("The product is: ")
        else:
            print(f"This produces {num_products} products:")
            self.display_mol_tuple(products)
            print(f"\nPick product #{choice + 1}:")

    # ---------------------------------------------------------

    # -------------------- For playground mode ----------------
    def display_2d_mol_tuple(self, products: Mol2dTuple) -> None:
        for index, scenario in enumerate(products, start=1):
            print(f"Scenario #{index}")
            self.display_mol_tuple(scenario)
            print("\n")

    def display_compatible_reactions(
        self, start_mol: Mol, possible_reactions: list[Reaction]
    ) -> None:
        self.display_mol(start_mol)
        print(f"Smiles of molecule: {Chem.MolToSmiles(start_mol)}")
        if possible_reactions:
            print("\nCompatible reactions are listed below. Choose a reaction:\n")
            for i, reaction in enumerate(possible_reactions, start=1):
                print(f"[{i}] {reaction.name}")
                print(reaction.description)
                print("\n")
        else:
            print("\nThere are no compatible reactions.\n")

    def prompt_for_scenario_index(self, products: Mol2dTuple) -> int:
        choice = input(f"\nPick a scenario to analyze next (1 - {len(products)}): ")
        return int(choice) - 1

    def prompt_for_product_index(self, scenario: MolTuple) -> int:
        choice = input(f"\nPick a product to analyze next (1 - {len(scenario)}): ")
        return int(choice) - 1

    # ---------------------------------------------------------
