from rdkit.Chem.rdchem import Mol
from PIL.Image import Image as PILImage

MolTuple = tuple[Mol, ...]
Mol2dTuple = tuple[tuple[Mol, ...], ...]

SolverModeImageData = dict[int, list[PILImage]]
