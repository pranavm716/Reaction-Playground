from collections import deque
from typing import Deque

from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from website.config import MAX_NUM_SOLVER_STEPS
from website.datatypes import Mol2dTuple, MolTuple
from website.reaction import Reaction


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
    reactants: Mol | MolTuple, reaction: Reaction
) -> Mol2dTuple:
    """
    Performs a reaction on a molecule once, and returns a 2d tuple of products.
    """
    if isinstance(reactants, Mol):
        reactants = (reactants,)

    products: Mol2dTuple = tuple()
    for subreaction in reaction:
        products += subreaction.RunReactants(reactants)

    return generate_unique_products(products)


def generate_multi_step_product(start_mol: Mol, reaction: Reaction) -> MolTuple:
    """
    Recursively performs a reaction on a starting molecule and on its products until no more products can be formed.
    Returns the final tuple of products (always 1 x n).
    """
    products = _generate_multi_step_product(start_mol, reaction, tuple())
    return generate_unique_products((products,))[0]


def _generate_multi_step_product(
    start_mol: Mol, reaction: Reaction, products: MolTuple
) -> MolTuple:
    """
    Recursive sub-method for the method above.
    """
    single_step_product = generate_single_step_product(start_mol, reaction)
    if not single_step_product:
        return (start_mol,)

    for scenario in single_step_product:
        for p in scenario:
            products += _generate_multi_step_product(p, reaction, products)

    return products


def get_reactant_position_of_mol_in_reaction(mol: Mol, reaction: Reaction) -> int:
    """
    Returns the position of the mol in the reactants of any of the reaction's subreactions.
    Since all the subreactions have the same substructure pattern, we can return the first
    instance of when the mol is found.
    """
    for subreaction in reaction:
        for index, reactant in enumerate(subreaction.GetReactants()):
            if mol.HasSubstructMatch(reactant):
                return index
    raise ValueError(
        f"{reaction.name} is not valid for molecule with SMILES {Chem.MolToSmiles(mol)!r}."
    )


def find_possible_reactions(
    start_mol: Mol,
    all_reactions: list[Reaction],
    *,
    solver_mode: bool,
) -> list[Reaction]:
    """
    Returns the reactions that can be performed on the given molecule.
    """
    possible_reactions = []
    for reaction in all_reactions:
        if solver_mode and reaction.num_reactants > 1:
            continue  # Solver mode will omit any reactions with more than one reactant

        for subreaction in reaction:
            if subreaction.IsMoleculeReactant(start_mol):
                possible_reactions.append(reaction)
                break

    return possible_reactions


def find_synthetic_pathway(
    start_mol: Mol, target_mol: Mol, all_reactions: list[Reaction]
) -> tuple[bool, list[Reaction], list[int]]:
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
    # key = (reactant SMILES, product SMILES) -> value = (reaction, product_index) that converts
    # the reactant into the product. Needed to reconstruct the reaction pathway later.
    choices: dict[tuple[str, str], tuple[Reaction, int]] = {}

    # BFS loop
    dist = 0
    path_found = False
    while queue and dist < MAX_NUM_SOLVER_STEPS:
        for _ in range(len(queue)):
            cur = queue.popleft()
            if Chem.MolToSmiles(cur) == Chem.MolToSmiles(target_mol):
                path_found = True
                break

            possible_reactions = find_possible_reactions(
                start_mol=cur, all_reactions=all_reactions, solver_mode=True
            )
            all_reactions_products: list[MolTuple] = [
                generate_multi_step_product(cur, rxn) for rxn in possible_reactions
            ]

            # Loop through all compatible reactions
            for reaction, multi_step_product in zip(
                possible_reactions, all_reactions_products
            ):
                # Loop through all products of each compatible reaction
                for product_index, mol in enumerate(multi_step_product):
                    if Chem.MolToSmiles(mol) not in visited:
                        visited.add(Chem.MolToSmiles(mol))
                        queue.append(mol)

                        parent[Chem.MolToSmiles(mol)] = Chem.MolToSmiles(cur)
                        choices[(Chem.MolToSmiles(cur), Chem.MolToSmiles(mol))] = (
                            reaction,
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
