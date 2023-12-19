from unittest import mock
from unittest.mock import MagicMock

from rdkit import Chem

from program import Program


@mock.patch("builtins.input", side_effect=["1", "CCC=O", "b", "q"])
def test_playground_mode(mock_input, config):
    mock_cli = MagicMock()
    mock_cli.get_user_input.side_effect = mock_input
    program = Program(
        ui=mock_cli,
        start_mol=Chem.MolFromSmiles("CC[CH-]C1CC1"),
        target_mol=None,
        multi_step_react_mode=config.multi_step_react_mode,
        max_num_solver_steps=config.max_num_solver_steps,
        all_reactions_file_path=config.all_reactions_file_path,
        multiple_reactants_prompts=config.multiple_reactants_prompts,
    )
    program.run_program()

    assert mock_cli.get_user_input.call_count == 4
    display_compatible_reactions_args = (
        mock_cli.display_compatible_reactions.call_args_list
    )
    start_mols = [
        Chem.MolToSmiles(arg.args[0]) for arg in display_compatible_reactions_args
    ]
    assert start_mols == [
        "CC[CH-]C1CC1",
        "CCC(O)C(CC)C1CC1",
        "CC[CH-]C1CC1",
    ]


@mock.patch("builtins.input", side_effect=["n"])
def test_solver_mode(mock_input, all_reactions, config):
    mock_cli = MagicMock()
    mock_cli.get_user_input.side_effect = mock_input

    start_mol_smiles = "CCCO"
    target_mol_smiles = "CCC(=O)Cl"
    program = Program(
        ui=mock_cli,
        start_mol=Chem.MolFromSmiles(start_mol_smiles),
        target_mol=Chem.MolFromSmiles(target_mol_smiles),
        multi_step_react_mode=config.multi_step_react_mode,
        max_num_solver_steps=config.max_num_solver_steps,
        all_reactions_file_path=config.all_reactions_file_path,
        multiple_reactants_prompts=config.multiple_reactants_prompts,
    )
    program.run_program()

    assert mock_cli.display_solver_mode_intro.call_count == 1

    display_solver_step_args = mock_cli.display_solver_step.call_args_list
    current_mols = [Chem.MolToSmiles(arg.args[0]) for arg in display_solver_step_args]
    assert current_mols == [start_mol_smiles, "CCC(=O)O"]

    reactions = [arg.args[2] for arg in display_solver_step_args]
    assert reactions == [all_reactions[i] for i in [2, 11]]

    choices = [arg.args[4] for arg in display_solver_step_args]
    assert choices == [0, 0]

    # Target molecule
    assert Chem.MolToSmiles(mock_cli.display_mol.call_args.args[0]) == target_mol_smiles
