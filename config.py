from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol
from reaction import Reaction

MolTuple = tuple[Mol, ...]
Mol2dTuple = tuple[tuple[Mol, ...], ...]

MULTI_STEP_REACT_MODE = True
MAX_SOLVER_STEPS = 15


def read_all_reactions_from_file() -> list[Reaction]:
    """
    Returns a list of all available reactions.
    """
    all_reactions_file_path = "./all_reactions.txt"
    with open(all_reactions_file_path, "r") as f:
        lines = [l.strip() for l in f.readlines()]

    all_reactions = []
    i = 0
    while i < len(lines):
        reaction_name = lines[i]
        i += 1
        reaction_list = []
        while lines[i].startswith("^"):
            reaction_list.append(AllChem.ReactionFromSmarts(lines[i][1:]))
            i += 1
        reactants_side = lines[i - 1].index(">")
        num_reactants = lines[i - 1][1:reactants_side].count(".") + 1
        reaction_description = lines[i]
        i += 2

        all_reactions.append(
            Reaction(
                reaction_name,
                reaction_list,
                reaction_description,
                num_reactants,
            )
        )

    return all_reactions
