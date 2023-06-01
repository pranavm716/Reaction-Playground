from config import Mol2dTuple
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from ui import UI


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

    def display_2d_mol_tuple(self, products: Mol2dTuple) -> None:
        for index, scenario in enumerate(products, start=1):
            print(f"Scenario #{index}")
            display(Draw.MolsToGridImage(list(scenario), molsPerRow=len(scenario)))
            print("\n")
