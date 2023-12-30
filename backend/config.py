import pathlib

# Backend data file paths
ALL_REACTIONS_FILE_PATH = (
    pathlib.Path(__file__).parent.resolve() / "data" / "all_reactions.yaml"
)
ALL_SUBSTRUCTURES_FILE_PATH = (
    pathlib.Path(__file__).parent.resolve() / "data" / "substructures.yaml"
)

# Frontend data file paths
STATIC_DIR = pathlib.Path(__file__).parent.resolve() / "static"
TEMPLATES_DIR = pathlib.Path(__file__).parent.resolve() / "templates"

# Computation configs
MULTI_STEP_REACT_MODE = True
MAX_NUM_SOLVER_STEPS = 15