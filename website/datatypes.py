from rdkit.Chem.rdchem import Mol
from PIL.Image import Image as PILImage

MolTuple = tuple[Mol, ...]
Mol2dTuple = tuple[tuple[Mol, ...], ...]

SmilesTuple = tuple[str, ...]
Smiles2dTuple = tuple[tuple[str, ...], ...]

# Mapping of the substructure name to its substructure patterns
SubstructureDict = dict[str, list[Mol]]

# Mapping of the step number to a list of images of the products generated
# by running the corresponding reaction for that step
SolverModeImageData = dict[int, list[PILImage]]
