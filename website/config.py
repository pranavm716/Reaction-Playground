import pathlib

MULTI_STEP_REACT_MODE = True
MAX_NUM_SOLVER_STEPS = 15
ALL_REACTIONS_FILE_PATH = (
    pathlib.Path(__file__).parent.resolve() / "data" / "all_reactions.yaml"
)
ALL_SUBSTRUCTURES_FILE_PATH = (
    pathlib.Path(__file__).parent.resolve() / "data" / "substructures.yaml"
)
