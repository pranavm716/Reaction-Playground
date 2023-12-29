import io
import pathlib
from unittest.mock import patch

import pytest
from rdkit import Chem, rdBase

from website.computations import (
    copy_mol,
    find_possible_reactions,
    find_synthetic_pathway,
    generate_multi_step_product,
    generate_single_step_product,
    generate_unique_products,
    get_reactant_position_of_mol_in_reaction,
)
from website.config import ALL_REACTIONS_FILE_PATH
from website.datatypes import Mol2dTuple, SmilesTuple, Smiles2dTuple
from website.reaction import Reaction, read_all_reactions_from_file


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
            "CO",
            1,  # DMP oxidation
            (("C=O",),),
        ],
        [
            "CCCO",
            5,  # PBr3 bromination of alcohols
            (("CCCBr",),),
        ],
        [
            r"C/C=C/C=C(\CC)C1CC1",
            6,  # Ozonolysis
            (("CCC(=O)C1CC1", "C/C=C/C=O"), ("CC=O", r"CC/C(=C\C=O)C1CC1")),
        ],
        [
            ("CC[CH-]C1CC1", "CCC=O"),
            13,  # Grignard reaction
            (("CCC(O)C(CC)C1CC1",),),
        ],
        [
            ("CCC[C-]", "N#CC1CCCCC1"),
            13,  # Grignard reaction
            (("CCCCC(=O)C1CCCCC1",),),
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
    else:
        raise TypeError("reactants_smiles must be of type str or tuple[str].")

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
        [
            "C=C(C#CCC1CC1)/C=C/C(C#CC/C=C/C)=C(/C)CC",
            6,
            ("C=O", "CCC(C)=O", "CC=O", "O=CC(=O)C(=O)O", "O=CCC(=O)O", "O=C(O)CC1CC1"),
        ],  # Ozonolysis,
        [
            "C=C(C#CCC1CC1)/C=C/C(C#CC/C=C/C)=C(/C)CC",
            8,
            (
                "CCC(C)(O)C(C#CCCC(C)O)C(O)CC(C)(O)C#CCC1CC1",
                "CCC(C)(O)C(C#CCCC(C)O)CC(O)C(C)(O)C#CCC1CC1",
                "CCC(C)C(O)(C#CCCC(C)O)C(O)CC(C)(O)C#CCC1CC1",
                "CCC(C)C(O)(C#CCCC(C)O)CC(O)C(C)(O)C#CCC1CC1",
                "CCC(O)CC#CC(C(O)CC(C)(O)C#CCC1CC1)C(C)(O)CC",
                "CCC(O)CC#CC(CC(O)C(C)(O)C#CCC1CC1)C(C)(O)CC",
                "CCC(O)CC#CC(O)(C(C)CC)C(O)CC(C)(O)C#CCC1CC1",
                "CCC(O)CC#CC(O)(CC(O)C(C)(O)C#CCC1CC1)C(C)CC",
            ),
        ],  # OM/DM
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
        ["C#CCC(O)C#C", 13, None],  # Grignard reaction
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


@pytest.mark.parametrize(
    ["start_mol_smiles", "solver_mode", "possible_reactions_indices"],
    [
        [r"C/C(C)=C\C=O", True, [2, 3, 4, 6, 8, 9]],
        [r"C/C(C)=C\C=O", False, [2, 3, 4, 6, 8, 9, 13, -1]],
        ["COC(=O)C1C#CCC1", False, [0, 4, 6]],
    ],
)
def test_find_possible_reactions(
    start_mol_smiles: str,
    solver_mode: bool,
    possible_reactions_indices: list[int],
    all_reactions: list[Reaction],
):
    start_mol = Chem.MolFromSmiles(start_mol_smiles)
    possible_reactions = find_possible_reactions(
        start_mol, all_reactions, solver_mode=solver_mode
    )
    correct_possible_reactions = [all_reactions[i] for i in possible_reactions_indices]
    assert possible_reactions == correct_possible_reactions


@pytest.mark.parametrize(
    [
        "start_mol_smiles",
        "target_mol_smiles",
        "path_found",
        "reaction_pathway_indices",
        "choice_pathway",
    ],
    [
        [
            "CC(=O)OC1=CCC(=O)C1",
            "OCC1CC=C(O)C1",
            True,
            [3, 5, 7, 0, 4],
            [0, 0, 0, 1, 0],
        ],
        [
            "OC1CC1C2C#CC2",
            "O=CCC1CC1=O",
            True,
            [6, 4, 1],
            [1, 0, 0],
        ],
        [
            "OC1CC1C2C#CC2",
            "O=CCC1CC1O",
            False,
            [],
            [],
        ],
        [
            "CCCC",
            "CCC#N",
            False,
            [],
            [],
        ],
        pytest.param(
            "C=C(C#CCC1CC1)/C=C/C(C#CC/C=C/C)=C(/C)CC",
            "CCC(C)C(O)(CC(=O)Cl)C(CC(C)(O)CC(=O)Cl)C(=O)Cl",
            True,
            [8, 6, 4, 5, 7, 0, 11],
            [2, 1, 0, 0, 0, 0, 0],
            marks=pytest.mark.slow,
        ),
    ],
)
def test_find_synthetic_pathway(
    start_mol_smiles: str,
    target_mol_smiles: str,
    path_found: bool,
    reaction_pathway_indices: list[int],
    choice_pathway: list[int],
    all_reactions: list[Reaction],
):
    start_mol = Chem.MolFromSmiles(start_mol_smiles)
    target_mol = Chem.MolFromSmiles(target_mol_smiles)

    (
        retrieved_path_found,
        retrieved_reaction_pathway,
        retrieved_choice_pathway,
    ) = find_synthetic_pathway(start_mol, target_mol, all_reactions)

    assert retrieved_path_found == path_found
    assert retrieved_reaction_pathway == [
        all_reactions[i] for i in reaction_pathway_indices
    ]
    assert retrieved_choice_pathway == choice_pathway


def test_all_reactions_are_configured_properly():
    # Ensure that rdkit does not throw any warnings
    rdBase.LogToPythonStderr()
    with patch("sys.stderr", new_callable=io.StringIO) as s:
        all_reactions_file_path = (
            pathlib.Path(__file__).resolve().parent.parent
            / "website"
            / ALL_REACTIONS_FILE_PATH
        )
        read_all_reactions_from_file(all_reactions_file_path)

    assert not s.getvalue()
