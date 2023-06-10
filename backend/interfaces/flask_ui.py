from flask import flash, request, session
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdchem import Mol

from backend.datatypes import Mol2dTuple, MolTuple
from backend.reaction import Reaction


class FlaskUI:
    """A concrete implementation of the UI class. This class contains methods
    to parse user input from the flask frontend and send image data to the flask
    frontend.
    """

    def save_mol(self, mol: Mol, file_path: str) -> None:
        img = Draw.MolToImage(mol)
        img.save(file_path, "PNG")

    def save_mol_tuple(self, mols: MolTuple, file_paths: list[str]) -> None:
        for mol, file_path in zip(mols, file_paths):
            self.save_mol(mol, file_path)

    # def get_user_input(self, prompt: str = "") -> str:
    #     return input(prompt)

    # def print(self, text: str = "") -> None:
    #     """
    #     Printing is taken care of in the frontend HTML.
    #     """
    #     pass

    # # -------------------- For playground mode ----------------
    # def display_2d_mol_tuple(self, products: Mol2dTuple) -> None:
    #     for index, scenario in enumerate(products, start=1):
    #         print(f"Scenario #{index}")
    #         self.save_mol_tuple(scenario)
    #         print("\n")

    # def display_compatible_reactions(
    #     self, start_mol: Mol, possible_reactions: list[Reaction]
    # ) -> None:
    #     self.save_mol(start_mol)
    #     print(f"Smiles of molecule: {Chem.MolToSmiles(start_mol)}")
    #     if possible_reactions:
    #         print("\nCompatible reactions are listed below. Choose a reaction:\n")
    #         for i, reaction in enumerate(possible_reactions, start=1):
    #             print(f"[{i}] {reaction.name}")
    #             print(reaction.description)
    #             print("\n")
    #     else:
    #         print("\nThere are no compatible reactions.\n")

    # def prompt_for_scenario_index(self, products: Mol2dTuple) -> int:
    #     choice = input(f"\nPick a scenario to analyze next (1 - {len(products)}): ")
    #     return int(choice) - 1

    # def prompt_for_product_index(self, scenario: MolTuple) -> int:
    #     choice = input(f"\nPick a product to analyze next (1 - {len(scenario)}): ")
    #     return int(choice) - 1

    # # ---------------------------------------------------------
