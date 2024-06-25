import pytest
from rdkit import Chem

from tests.utils import (
    mol_2d_tuple_to_smiles,
    smiles_2d_tuple_to_mols,
    smiles_2d_tuples_match,
)
from backend.computations import (
    copy_mol,
    find_possible_reaction_keys,
    find_synthetic_pathway,
    generate_multi_step_product,
    generate_single_step_product,
    generate_unique_products,
    get_reactant_position_of_mol_in_reaction,
    get_substructure_classifications,
)
from backend.datatypes import SmilesTuple, Smiles2dTuple
from backend.mol_classification import MolClass
from backend.reaction import ReactionKey


class TestComputationUtils:
    def test_copy_mol(self):
        smiles = "c1ccc2c(c1)c3ccccc3[nH]2"
        mol = Chem.MolFromSmiles(smiles)
        copy = copy_mol(mol)
        copy_smiles = Chem.MolToSmiles(copy)
        assert Chem.CanonSmiles(smiles) == Chem.CanonSmiles(copy_smiles)

    def test_generate_unique_products(self):
        product_smiles = (
            ("C", "CCC", "C"),
            ("CCC", "CCC", "C"),
            ("N",),
            ("CCCO", "CCC#N"),
        )
        products = smiles_2d_tuple_to_mols(product_smiles)

        unique_products = generate_unique_products(products)
        unique_product_smiles = mol_2d_tuple_to_smiles(unique_products)

        assert smiles_2d_tuples_match(
            unique_product_smiles, (("N",), ("CCC#N", "CCCO"), ("C", "CCC"))
        )


class TestProductGeneration:
    @pytest.mark.parametrize(
        ["reactants_smiles", "reaction_key", "single_step_product_smiles"],
        [
            [
                "CCO",
                ReactionKey.h2cro4_oxidation,
                (("CC(=O)O",),),
            ],
            [
                "CO",
                ReactionKey.dmp_pcc_oxidation,
                (("C=O",),),
            ],
            [
                "CCCO",
                ReactionKey.pbr3_bromination,
                (("CCCBr",),),
            ],
            [
                r"C/C=C/C=C(\CC)C1CC1",
                ReactionKey.ozonolysis,
                (("CCC(=O)C1CC1", "C/C=C/C=O"), ("CC=O", r"CC/C(=C\C=O)C1CC1")),
            ],
            [
                ("CC[CH-]C1CC1", "CCC=O"),
                ReactionKey.grignard_reaction,
                (("CCC(O)C(CC)C1CC1",),),
            ],
            [
                ("CCC[C-]", "N#CC1CCCCC1"),
                ReactionKey.grignard_reaction,
                (("CCCCC(=O)C1CCCCC1",),),
            ],
        ],
    )
    def test_generate_single_step_product(
        self,
        reactants_smiles: str | SmilesTuple,
        reaction_key: ReactionKey,
        single_step_product_smiles: Smiles2dTuple,
    ):
        if isinstance(reactants_smiles, str):
            reactants = Chem.MolFromSmiles(reactants_smiles)
        elif isinstance(reactants_smiles, tuple):
            reactants = tuple(Chem.MolFromSmiles(s) for s in reactants_smiles)
        else:
            raise TypeError("reactants_smiles must be of type str or tuple[str].")

        single_step_product = generate_single_step_product(reactants, reaction_key)
        assert smiles_2d_tuples_match(
            single_step_product_smiles, mol_2d_tuple_to_smiles(single_step_product)
        )

    @pytest.mark.parametrize(
        ["reactant_smiles", "reaction_key", "multi_step_product_smiles"],
        [
            [
                "CC(=O)OC1=CCC(=O)C1",
                ReactionKey.lialh4_reduction,
                ("CCO", "OC1=CCC(O)C1"),
            ],
            ["C=CC", ReactionKey.om_dm, ("CC(C)O",)],
            ["C=CC", ReactionKey.hydroboration_oxidation, ("CCCO",)],
            [
                "C=C(C#CCC1CC1)/C=C/C(C#CC/C=C/C)=C(/C)CC",
                ReactionKey.ozonolysis,
                (
                    "C=O",
                    "CCC(C)=O",
                    "CC=O",
                    "O=CC(=O)C(=O)O",
                    "O=CCC(=O)O",
                    "O=C(O)CC1CC1",
                ),
            ],
            [
                "C=C(C#CCC1CC1)/C=C/C(C#CC/C=C/C)=C(/C)CC",
                ReactionKey.om_dm,
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
            ],
        ],
    )
    def test_generate_multi_step_product(
        self,
        reactant_smiles: str,
        reaction_key: ReactionKey,
        multi_step_product_smiles: SmilesTuple,
    ):
        start_mol = Chem.MolFromSmiles(reactant_smiles)
        multi_step_product = generate_multi_step_product(start_mol, reaction_key)
        assert set(multi_step_product_smiles) == set(
            Chem.MolToSmiles(p) for p in multi_step_product
        )


class TestReactionQueryMethods:
    @pytest.mark.parametrize(
        ["reactant_smiles", "reaction_key", "reactant_position"],
        [
            ["BrCC1CCCC1", ReactionKey.nacn_nitrile_synthesis, 0],
            ["C#CCC(O)C#C", ReactionKey.grignard_reaction, None],
        ],
    )
    def test_get_reactant_position(
        self,
        reactant_smiles: str,
        reaction_key: ReactionKey,
        reactant_position: int | None,
    ):
        mol = Chem.MolFromSmiles(reactant_smiles)

        if reactant_position is None:
            with pytest.raises(ValueError):
                get_reactant_position_of_mol_in_reaction(mol, reaction_key)
        else:
            assert (
                get_reactant_position_of_mol_in_reaction(mol, reaction_key)
                == reactant_position
            )

    @pytest.mark.parametrize(
        ["start_mol_smiles", "solver_mode", "expected_possible_reactions_keys"],
        [
            [
                r"C/C(C)=C\C=O",
                True,
                [
                    ReactionKey.h2cro4_oxidation,
                    ReactionKey.nabh4_reduction,
                    ReactionKey.lialh4_reduction,
                    ReactionKey.ozonolysis,
                    ReactionKey.om_dm,
                    ReactionKey.hydroboration_oxidation,
                ],
            ],
            [
                r"C/C(C)=C\C=O",
                False,
                [
                    ReactionKey.h2cro4_oxidation,
                    ReactionKey.nabh4_reduction,
                    ReactionKey.lialh4_reduction,
                    ReactionKey.ozonolysis,
                    ReactionKey.om_dm,
                    ReactionKey.hydroboration_oxidation,
                    ReactionKey.grignard_reaction,
                ],
            ],
            [
                "COC(=O)C1C#CCC1",
                False,
                [
                    ReactionKey.hydrolysis,
                    ReactionKey.lialh4_reduction,
                    ReactionKey.ozonolysis,
                ],
            ],
        ],
    )
    def test_find_possible_reaction_keys(
        self,
        start_mol_smiles: str,
        solver_mode: bool,
        expected_possible_reactions_keys: list[ReactionKey],
    ):
        start_mol = Chem.MolFromSmiles(start_mol_smiles)
        possible_reaction_keys = find_possible_reaction_keys(
            start_mol, solver_mode=solver_mode
        )
        assert possible_reaction_keys == expected_possible_reactions_keys


class TestFindSyntheticPathway:
    @pytest.mark.parametrize(
        [
            "start_mol_smiles",
            "target_mol_smiles",
            "path_found",
            "reaction_pathway_keys",
            "choice_pathway",
        ],
        [
            [
                "CC(=O)OC1=CCC(=O)C1",
                "OCC1CC=C(O)C1",
                True,
                [
                    ReactionKey.nabh4_reduction,
                    ReactionKey.pbr3_bromination,
                    ReactionKey.nacn_nitrile_synthesis,
                    ReactionKey.hydrolysis,
                    ReactionKey.lialh4_reduction,
                ],
                [0, 0, 0, 1, 0],
            ],
            [
                "OC1CC1C2C#CC2",
                "O=CCC1CC1=O",
                True,
                [
                    ReactionKey.ozonolysis,
                    ReactionKey.lialh4_reduction,
                    ReactionKey.dmp_pcc_oxidation,
                ],
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
                [
                    ReactionKey.om_dm,
                    ReactionKey.ozonolysis,
                    ReactionKey.lialh4_reduction,
                    ReactionKey.pbr3_bromination,
                    ReactionKey.nacn_nitrile_synthesis,
                    ReactionKey.hydrolysis,
                    ReactionKey.socl2_acid_chloride_synthesis,
                ],
                [2, 1, 0, 0, 0, 0, 0],
                marks=pytest.mark.slow,
            ),
        ],
    )
    def test_find_synthetic_pathway(
        self,
        start_mol_smiles: str,
        target_mol_smiles: str,
        path_found: bool,
        reaction_pathway_keys: list[ReactionKey],
        choice_pathway: list[int],
    ):
        start_mol = Chem.MolFromSmiles(start_mol_smiles)
        target_mol = Chem.MolFromSmiles(target_mol_smiles)

        (
            retrieved_path_found,
            retrieved_reaction_pathway_keys,
            retrieved_choice_pathway,
        ) = find_synthetic_pathway(start_mol, target_mol)

        assert retrieved_path_found == path_found
        assert retrieved_reaction_pathway_keys == reaction_pathway_keys
        assert retrieved_choice_pathway == choice_pathway


class TestSubstructureClassifications:
    @pytest.mark.parametrize(
        ["mol_smiles", "expected_substructures"],
        [
            ["N#CC1CCCCC1", [MolClass.nitrile]],
            ["CCCCC(=O)C1CCCCC1", [MolClass.ketone]],
            [
                "C#CCC(=O)O",
                [MolClass.carboxylic_acid, MolClass.terminal_alkyne],
            ],
            [
                "CC(C)(O)CCC(Br)C=O",
                [
                    MolClass.aldehyde,
                    MolClass.secondary_alkyl_bromide,
                    MolClass.tertiary_alcohol,
                ],
            ],
            [
                "CCCOCC=CC(=O)NC",
                [MolClass.internal_alkene, MolClass.ether, MolClass.secondary_amide],
            ],
            ["CCN(C)CC=O", [MolClass.tertiary_amine, MolClass.aldehyde]],
            [
                "CCOC(=O)CCC(=O)N(C)C",
                [MolClass.tertiary_amide, MolClass.ester],
            ],
            [
                "[C-]CC#CC(=O)Cl",
                [
                    MolClass.acid_chloride,
                    MolClass.internal_alkyne,
                    MolClass.carbon_nucleophile,
                ],
            ],
            [
                "NC(=O)CC(O)C(CBr)=C",
                [
                    MolClass.primary_alkyl_bromide,
                    MolClass.secondary_alcohol,
                    MolClass.primary_amide,
                    MolClass.terminal_alkene,
                ],
            ],
            [
                "NC(CO)CCNC1CCC1",
                [
                    MolClass.primary_alcohol,
                    MolClass.secondary_amine,
                    MolClass.primary_amine,
                ],
            ],
            ["C", []],
            # Single carbon molecules
            ["C=O", [MolClass.aldehyde]],
            ["CO", [MolClass.primary_alcohol]],
            ["C(=O)O", [MolClass.carboxylic_acid]],
            ["CBr", [MolClass.primary_alkyl_bromide]],
            ["C=C", [MolClass.terminal_alkene]],
            ["C#C", [MolClass.terminal_alkyne]],
        ],
    )
    def test_substructure_classifications(
        self,
        mol_smiles: str,
        expected_substructures: list[str],
    ):
        mol = Chem.MolFromSmiles(mol_smiles)
        assert set(get_substructure_classifications(mol)) == set(expected_substructures)
