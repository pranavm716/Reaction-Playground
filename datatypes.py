from rdkit.Chem.rdchem import Mol
from PIL.Image import Image as PILImage

MolTuple = tuple[Mol, ...]
Mol2dTuple = tuple[tuple[Mol, ...], ...]
SolverModeImageData = list[PILImage | tuple[list[PILImage], int]]
