from pydantic import BaseModel


class MolImageMetadata(BaseModel):
    smiles: str
    encoding: str


# Mapping of the step number to a list of smiles + image encodings of the products generated
# by running the corresponding reaction for that step
SolverModeImageData = dict[int, list[MolImageMetadata]]


class SolverModeResponse(BaseModel):
    path_found: bool
    num_steps: int
    reaction_names: list[str]
    choice_pathway: list[int]
    starting_encoding: str
    target_encoding: str
    solver_image_metadata: SolverModeImageData
