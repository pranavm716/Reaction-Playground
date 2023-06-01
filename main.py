from config import read_config
from program import Program
from ui import GoogleColabUI


def main() -> None:
    # Read the configuration settings from a json file
    config = read_config("./config.json")

    # Prompt the user for the start and target molecules
    google_colab_ui = GoogleColabUI()
    start_mol, target_mol = google_colab_ui.prompt_for_start_and_target_molecule()

    program = Program(
        ui=google_colab_ui,
        start_mol=start_mol,
        target_mol=target_mol,
        multi_step_react_mode=config.multi_step_react_mode,
        max_num_solver_steps=config.max_num_solver_steps,
        all_reactions_file_path=config.all_reactions_file_path,
    )
    program.run()


if __name__ == "__main__":
    main()
