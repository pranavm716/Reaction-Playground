from rdkit import Chem

from website.datatypes import Mol2dTuple, Smiles2dTuple


def mol_2d_tuple_to_smiles(mol_2d_tuple: Mol2dTuple) -> Smiles2dTuple:
    return tuple(
        tuple(Chem.MolToSmiles(s) for s in scenarios) for scenarios in mol_2d_tuple
    )


def smiles_2d_tuple_to_mols(smiles_2d_tuple: Smiles2dTuple) -> Mol2dTuple:
    return tuple(
        tuple(Chem.MolFromSmiles(s) for s in scenarios) for scenarios in smiles_2d_tuple
    )


def smiles_2d_tuples_match(tuple1: Smiles2dTuple, tuple2: Smiles2dTuple) -> bool:
    """
    Helper method to check if two smiles 2d tuples are the same. Takes into consideration that
    the order of the inner and outer tuples does not matter.
    """
    if len(tuple1) != len(tuple2):
        return False

    set1, set2 = set(tuple1), set(tuple2)
    for inner1 in set1:
        inner_set1 = set(inner1)
        for inner2 in set2:
            inner_set2 = set(inner2)
            if inner_set1 == inner_set2:
                break
        else:
            return False
    return True
