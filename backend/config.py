import pathlib
import os

# Backend data file paths
ALL_REACTIONS_FILE_PATH = (
    pathlib.Path(__file__).parent.resolve() / "data" / "all_reactions.yaml"
)
ALL_SUBSTRUCTURES_FILE_PATH = (
    pathlib.Path(__file__).parent.resolve() / "data" / "substructures.yaml"
)

# Computation configs
MULTI_STEP_REACT_MODE = True
MAX_NUM_SOLVER_STEPS = 15

# Frontend
REACT_JS_REQUEST_ORIGIN = os.getenv("FRONTEND_URL", "http://localhost:3000")
