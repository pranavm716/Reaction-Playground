from config import read_config
from program import Program
from ui import CLI


def main() -> None:
    # Read the configuration settings from a json file
    config = read_config("./config.json")

    if config.disable_rdkit_warnings:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.warning")

    # Prompt the user for the start and target molecules
    cli = CLI()
    start_mol, target_mol = cli.prompt_for_start_and_target_molecule()

    program = Program(
        ui=cli,
        start_mol=start_mol,
        target_mol=target_mol,
        multi_step_react_mode=config.multi_step_react_mode,
        max_num_solver_steps=config.max_num_solver_steps,
        all_reactions_file_path=config.all_reactions_file_path,
        multiple_reactants_prompts=config.multiple_reactants_prompts,
    )
    program.run_program()


if __name__ == "__main__":
    main()
