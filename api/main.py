from fastapi import FastAPI
from rdkit import Chem
from starlette.middleware.cors import CORSMiddleware

from backend.computations import (
    ALL_REACTIONS,
    get_substructure_classifications,
)
from backend.config import REACT_JS_REQUEST_ORIGIN
from backend.reaction import ReactionKey, ReactionDict
from api.playground_router import router as playground_router

app = FastAPI()

# Add middleware to interface with React frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=[REACT_JS_REQUEST_ORIGIN],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/reactions")
def all_reactions() -> ReactionDict:
    """Returns a list of all available reactions."""
    return ALL_REACTIONS


@app.get("/reactions/{reaction_key}", response_model=ReactionDict)
def get_reaction(reaction_key: ReactionKey) -> ReactionDict:
    """Gets an individual reaction by its reaction key."""
    return {reaction_key: ALL_REACTIONS[reaction_key]}


@app.get("/mol/classifications")
def get_mol_classifications(mol_smiles: str) -> list[str]:
    return get_substructure_classifications(Chem.MolFromSmiles(mol_smiles))


app.include_router(playground_router)