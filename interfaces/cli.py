from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdchem import Mol

from datatypes import Mol2dTuple, MolTuple
from reaction import Reaction
from interfaces.ui import UI


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
        # img.show()
        import os

        fp = os.path.join(
            os.path.dirname(__file__),
            "../frontend/website/static/images/dynamic/intro.png",
        )
        img.save(fp, "PNG")

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
