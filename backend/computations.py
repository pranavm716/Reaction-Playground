from collections import deque
from typing import Deque, Generator

from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from backend.config import (
    MAX_NUM_SOLVER_STEPS,
    ALL_REACTIONS_FILE_PATH,
    ALL_SUBSTRUCTURES_FILE_PATH,
)
from backend.datatypes import Mol2dTuple, MolTuple
from backend.mol_classification import (
    SubstructureDict,
    read_all_substructures_from_file,
)
from backend.reaction import (
    ReactionKey,
    ReactionDict,
    read_all_reactions_from_file,
    Reaction,
)

ALL_REACTIONS: ReactionDict = read_all_reactions_from_file(ALL_REACTIONS_FILE_PATH)
ALL_SUBSTRUCTURES: SubstructureDict = read_all_substructures_from_file(
    ALL_SUBSTRUCTURES_FILE_PATH
)


def copy_mol(mol: Mol) -> Mol:
    return Chem.MolFromSmiles(Chem.MolToSmiles(mol))


def generate_unique_products(products: Mol2dTuple) -> Mol2dTuple:
    """
    Returns only the unique products in a 2d mol tuple.
    """

    # Convert to SMILES
    smiles = tuple(
        tuple(Chem.MolToSmiles(s) for s in scenarios) for scenarios in products
    )
    inner_unique = tuple(tuple(sorted(set(p))) for p in smiles)
    outer_unique = tuple(sorted(set(inner_unique)))

    # Convert back to Mols
    return tuple(
        tuple(Chem.MolFromSmiles(s) for s in scenarios) for scenarios in outer_unique
    )


def generate_single_step_product(
    reactants: Mol | MolTuple, reaction_key: ReactionKey
) -> Mol2dTuple:
    """
    Performs a reaction on a molecule once, and returns a 2d tuple of products.
    """
    if isinstance(reactants, Mol):
        reactants = (reactants,)

    products: Mol2dTuple = tuple()
    for subreaction in ALL_REACTIONS[reaction_key]:
        products += subreaction.RunReactants(reactants)

    return generate_unique_products(products)


def generate_multi_step_product(start_mol: Mol, reaction_key: ReactionKey) -> MolTuple:
    """
    Recursively performs a reaction on a starting molecule and on its products until no more products can be formed.
    Returns the final tuple of products (always 1 x n).
    """

    # Use a set to keep track of which products have already been processed
    processed: set[str] = set()

    def _generate_multi_step_product(current_mol: Mol) -> Generator[Mol, None, None]:
        """
        Recursive sub-method for the method above.
        """
        if (smiles := Chem.MolToSmiles(current_mol)) in processed:
            return
        processed.add(smiles)

        single_step_product = generate_single_step_product(current_mol, reaction_key)
        if not single_step_product:
            yield current_mol
            return

        for scenario in single_step_product:
            for p in scenario:
                yield from _generate_multi_step_product(p)

    products = _generate_multi_step_product(start_mol)
    return tuple(sorted(set(products), key=lambda x: Chem.MolToSmiles(x)))


def get_reactant_position_of_mol_in_reaction(
    mol: Mol, reaction_key: ReactionKey
) -> int:
    """
    Returns the position of the mol in the reactants of any of the reaction's subreactions.
    Since all the subreactions have the same substructure pattern, we can return the first
    instance of when the mol is found.
    """
    reaction = ALL_REACTIONS[reaction_key]
    for subreaction in reaction:
        for index, reactant in enumerate(subreaction.GetReactants()):
            if mol.HasSubstructMatch(reactant):
                return index
    raise ValueError(
        f"{reaction.name} is not valid for molecule with SMILES {Chem.MolToSmiles(mol)!r}."
    )


def find_possible_reaction_keys(
    start_mol: Mol, *, solver_mode: bool
) -> list[ReactionKey]:
    """
    Returns the keys of the reactions that can be performed on the given molecule.
    """

    possible_reaction_keys = []
    for key, reaction in ALL_REACTIONS.items():
        if solver_mode and reaction.num_reactants > 1:
            continue  # Solver mode will omit any reactions with more than one reactant

        for subreaction in reaction:
            if subreaction.IsMoleculeReactant(start_mol):
                possible_reaction_keys.append(key)
                break

    return possible_reaction_keys


def get_reactions_from_keys(reaction_keys: list[ReactionKey]) -> list[Reaction]:
    return [ALL_REACTIONS[key] for key in reaction_keys]


def get_substructure_classifications(mol: Mol) -> list[str]:
    """
    Returns a list of classifications (essentially functional groups) that the given mol falls under.
    """
    return [
        name
        for name, patterns in ALL_SUBSTRUCTURES.items()
        if any(mol.HasSubstructMatch(pattern) for pattern in patterns)
    ]


def find_synthetic_pathway(
    start_mol: Mol, target_mol: Mol
) -> tuple[bool, list[ReactionKey], list[int]]:
    """
    Auto-solver, utilizes multi step products.
    Performs breadth first search (BFS) to find the SHORTEST possible reaction pathway to the target molecule.
    """
    if not target_mol:
        raise ValueError("Target molecule must be present to use the solver mode.")

    # Set of visited molecule SMILES - used to prevent infinite loops
    visited: set[str] = set()
    queue: Deque[Mol] = deque()

    # Add the start_node to the queue and visited list
    queue.append(start_mol)
    visited.add(Chem.MolToSmiles(start_mol))

    # The format of this dictionary is:
    # key = product SMILES -> value = reactant SMILES
    parent: dict[str, str | None] = {
        Chem.MolToSmiles(start_mol): None,  # startMol has no previous reactants
    }

    # The format of this dictionary is:
    # key = (reactant SMILES, product SMILES) -> value = (reaction key, product_index) that converts
    # the reactant into the product. Needed to reconstruct the reaction pathway later.
    choices: dict[tuple[str, str], tuple[ReactionKey, int]] = {}

    # BFS loop
    dist = 0
    path_found = False
    while queue and dist < MAX_NUM_SOLVER_STEPS:
        for _ in range(len(queue)):
            cur = queue.popleft()
            if Chem.MolToSmiles(cur) == Chem.MolToSmiles(target_mol):
                path_found = True
                break

            possible_reaction_keys = find_possible_reaction_keys(
                start_mol=cur, solver_mode=True
            )
            all_reactions_products: list[MolTuple] = [
                generate_multi_step_product(cur, key) for key in possible_reaction_keys
            ]

            # Loop through all compatible reactions
            for reaction_key, multi_step_product in zip(
                possible_reaction_keys, all_reactions_products
            ):
                # Loop through all products of each compatible reaction
                for product_index, mol in enumerate(multi_step_product):
                    if Chem.MolToSmiles(mol) not in visited:
                        visited.add(Chem.MolToSmiles(mol))
                        queue.append(mol)

                        parent[Chem.MolToSmiles(mol)] = Chem.MolToSmiles(cur)
                        choices[(Chem.MolToSmiles(cur), Chem.MolToSmiles(mol))] = (
                            reaction_key,
                            product_index,
                        )
        if path_found:
            break
        dist += 1

    if not path_found:
        return False, [], []

    # Reconstructing the reaction pathway
    reaction_pathway = []
    choice_pathway = []
    target_smiles = Chem.MolToSmiles(target_mol)

    while parent[target_smiles] is not None:
        choice = choices[(parent[target_smiles], target_smiles)]
        reaction, product_index = choice
        reaction_pathway.append(reaction)
        choice_pathway.append(product_index)
        target_smiles = parent[target_smiles]

    reaction_pathway.reverse()
    choice_pathway.reverse()

    return path_found, reaction_pathway, choice_pathway
