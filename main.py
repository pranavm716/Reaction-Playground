from config import read_all_reactions_from_file
from ui import GoogleColabUI
from program import Program


def main() -> None:
    google_colab_ui = GoogleColabUI()
    start_mol, target_mol = google_colab_ui.prompt_for_start_and_target_molecule()
    program = Program(
        ui=google_colab_ui,
        all_reactions=read_all_reactions_from_file(),
        start_mol=start_mol,
        target_mol=target_mol,
    )
    program.run()


if __name__ == "__main__":
    main()
