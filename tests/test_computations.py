from rdkit import Chem

from computations import (
    copy_mol,
    generate_single_step_product,
    generate_unique_products,
)
from datatypes import Mol2dTuple


def mol_2d_tuple_to_smiles(mol_2d_tuple: Mol2dTuple) -> tuple[tuple[str, ...], ...]:
    return tuple(
        tuple(Chem.MolToSmiles(s) for s in scenarios) for scenarios in mol_2d_tuple
    )


def smiles_2d_tuple_to_mols(smiles_2d_tuple: tuple[tuple[str, ...], ...]) -> Mol2dTuple:
    return tuple(
        tuple(Chem.MolFromSmiles(s) for s in scenarios) for scenarios in smiles_2d_tuple
    )


def test_copy_mol():
    mol = Chem.MolFromSmiles("CCC(C#N)C1CC1")
    copy = copy_mol(mol)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(copy)


def test_generate_unique_products():
    product_smiles = (("C", "CCC", "C"), ("CCC", "CCC", "C"), ("N",), ("CCCO", "CCC#N"))
    products = smiles_2d_tuple_to_mols(product_smiles)

    unique_products = generate_unique_products(products)
    unique_product_smiles = mol_2d_tuple_to_smiles(unique_products)

    assert set(unique_product_smiles) == set((("N",), ("CCC#N", "CCCO"), ("C", "CCC")))


def test_generate_single_step_product(all_reactions):
    pass
