from reaction import read_all_reactions_from_file
from website import create_app
import pathlib

MULTI_STEP_REACT_MODE = True
MAX_NUM_SOLVER_STEPS = 15
DISABLE_RDKIT_WARNINGS = True

ALL_REACTIONS_FILE_PATH = pathlib.Path("./all_reactions.yaml")
ALL_REACTIONS = read_all_reactions_from_file(ALL_REACTIONS_FILE_PATH)


def main():
    if DISABLE_RDKIT_WARNINGS:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.warning")

    app = create_app()
    app.run(debug=True)


if __name__ == "__main__":
    main()
