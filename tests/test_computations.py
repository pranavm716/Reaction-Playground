import pytest
from rdkit import Chem

from computations import (
    copy_mol,
    generate_multi_step_product,
    generate_single_step_product,
    generate_unique_products,
    get_reactant_position_of_mol_in_reaction,
)
from datatypes import Mol2dTuple
from reaction import Reaction

SmilesTuple = tuple[str, ...]
Smiles2dTuple = tuple[tuple[str, ...], ...]


def mol_2d_tuple_to_smiles(mol_2d_tuple: Mol2dTuple) -> Smiles2dTuple:
    return tuple(
        tuple(Chem.MolToSmiles(s) for s in scenarios) for scenarios in mol_2d_tuple
    )


def smiles_2d_tuple_to_mols(smiles_2d_tuple: Smiles2dTuple) -> Mol2dTuple:
    return tuple(
        tuple(Chem.MolFromSmiles(s) for s in scenarios) for scenarios in smiles_2d_tuple
    )


def smiles_2d_tuples_match(tuple1: Smiles2dTuple, tuple2: Smiles2dTuple) -> bool:
    """
    Helper method to check if two smiles 2d tuples are the same. Takes into consideration that
    the order of the inner and outer tuples does not matter.
    """
    if len(tuple1) != len(tuple2):
        return False

    set1, set2 = set(tuple1), set(tuple2)
    for inner1 in set1:
        inner_set1 = set(inner1)
        for inner2 in set2:
            inner_set2 = set(inner2)
            if inner_set1 == inner_set2:
                break
        else:
            return False
    return True


def test_copy_mol():
    smiles = "c1ccc2c(c1)c3ccccc3[nH]2"
    mol = Chem.MolFromSmiles(smiles)
    copy = copy_mol(mol)
    copy_smiles = Chem.MolToSmiles(copy)
    assert Chem.CanonSmiles(smiles) == Chem.CanonSmiles(copy_smiles)


def test_generate_unique_products():
    product_smiles = (("C", "CCC", "C"), ("CCC", "CCC", "C"), ("N",), ("CCCO", "CCC#N"))
    products = smiles_2d_tuple_to_mols(product_smiles)

    unique_products = generate_unique_products(products)
    unique_product_smiles = mol_2d_tuple_to_smiles(unique_products)

    assert smiles_2d_tuples_match(
        unique_product_smiles, (("N",), ("CCC#N", "CCCO"), ("C", "CCC"))
    )


@pytest.mark.parametrize(
    ["reactants_smiles", "reaction_index", "single_step_product_smiles"],
    [
        [
            "CCO",
            2,  # H2CrO4 oxidation
            (("CC(=O)O",),),
        ],
        [
            r"C/C=C/C=C(\CC)C1CC1",
            6,  # Ozonolysis
            (("CCC(=O)C1CC1", "C/C=C/C=O"), ("CC=O", r"CC/C(=C\C=O)C1CC1")),
        ],
        [
            ("CC[CH-]C1CC1", "CCC=O"),
            12,  # Grignard reaction
            (("CCC(O)C(CC)C1CC1",),),
        ],
    ],
)
def test_generate_single_step_product(
    reactants_smiles: str | SmilesTuple,
    reaction_index: int,
    single_step_product_smiles: Smiles2dTuple,
    all_reactions: list[Reaction],
):
    reaction = all_reactions[reaction_index]
    if isinstance(reactants_smiles, str):
        reactants = Chem.MolFromSmiles(reactants_smiles)
    elif isinstance(reactants_smiles, tuple):
        reactants = tuple(Chem.MolFromSmiles(s) for s in reactants_smiles)
    single_step_product = generate_single_step_product(reactants, reaction)
    assert smiles_2d_tuples_match(
        single_step_product_smiles, mol_2d_tuple_to_smiles(single_step_product)
    )


@pytest.mark.parametrize(
    ["reactant_smiles", "reaction_index", "multi_step_product_smiles"],
    [
        ["CC(=O)OC1=CCC(=O)C1", 4, ("CCO", "OC1=CCC(O)C1")],  # LiAlH4 reduction
        ["C=CC", 8, ("CC(C)O",)],  # OM/DM
        ["C=CC", 9, ("CCCO",)],  # BH3/[O]
    ],
)
def test_generate_multi_step_product(
    reactant_smiles: str,
    reaction_index: int,
    multi_step_product_smiles: SmilesTuple,
    all_reactions: list[Reaction],
):
    start_mol = Chem.MolFromSmiles(reactant_smiles)
    reaction = all_reactions[reaction_index]
    multi_step_product = generate_multi_step_product(start_mol, reaction)
    assert set(multi_step_product_smiles) == set(
        Chem.MolToSmiles(p) for p in multi_step_product
    )


@pytest.mark.parametrize(
    ["reactant_smiles", "reaction_index", "reactant_position"],
    [
        ["BrCC1CCCC1", 7, 0],  # NaCN nitrile synthesis
        ["CCCO", -1, 2],  # A made up reaction
        ["C#CCC(O)C#C", 12, None],  # Grignard reaction
    ],
)
def test_get_reactant_position(
    reactant_smiles: str,
    reaction_index: int,
    reactant_position: int | None,
    all_reactions: list[Reaction],
):
    mol = Chem.MolFromSmiles(reactant_smiles)
    reaction = all_reactions[reaction_index]

    if reactant_position is None:
        with pytest.raises(ValueError):
            get_reactant_position_of_mol_in_reaction(mol, reaction)
    else:
        assert (
            get_reactant_position_of_mol_in_reaction(mol, reaction) == reactant_position
        )
