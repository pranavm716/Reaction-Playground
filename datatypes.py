from rdkit.Chem.rdchem import Mol
from PIL.Image import Image as PILImage

MolTuple = tuple[Mol, ...]
Mol2dTuple = tuple[tuple[Mol, ...], ...]

MolImageData = tuple[PILImage, str]
SolverModeImageData = list[tuple[list[MolImageData], int]]
