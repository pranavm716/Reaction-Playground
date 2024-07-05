from rdkit.Chem.rdchem import Mol

MolTuple = tuple[Mol, ...]
Mol2dTuple = tuple[tuple[Mol, ...], ...]

SmilesTuple = tuple[str, ...]
Smiles2dTuple = tuple[tuple[str, ...], ...]
