from collections import deque

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol

from datatypes import Mol2dTuple
from datatypes import MolTuple
from reaction import Reaction


def generate_unique_products(products: Mol2dTuple) -> Mol2dTuple:
    """
    Returns only the unique products in a 2d mol tuple.
    """

    # Convert to SMILES
    smiles = tuple(
        tuple(Chem.MolToSmiles(s) for s in scenarios) for scenarios in products
    )
    inner_unique = tuple(set(p) for p in smiles)
    outer_unique = tuple(set(inner_unique))

    # Convert back to Mols
    return tuple(
        tuple(Chem.MolFromSmiles(s) for s in scenarios) for scenarios in outer_unique
    )


def generate_single_step_product(
    start_mol: Mol,
    reaction: Reaction,
    additional_reactants: MolTuple | None = None,
) -> Mol2dTuple:
    """
    Performs a reaction on a molecule once, and returns a 2d tuple of products.
    """
    all_reactants: MolTuple = (start_mol,)
    if additional_reactants:
        all_reactants += additional_reactants

    products: Mol2dTuple = tuple()
    for r in reaction.reactions_list:
        products += r.RunReactants(all_reactants)

    return generate_unique_products(products)


def generate_multi_step_product(start_mol: Mol, reaction: Reaction) -> MolTuple:
    """
    Recursively performs a reaction on a starting molecule and on its products until no more products can be formed.
    Returns the final tuple of products (always 1 x n).
    """
    products = _generate_multi_step_product_recursion(start_mol, reaction, tuple())
    return generate_unique_products((products,))[0]


def _generate_multi_step_product_recursion(
    start_mol: Mol, reaction: Reaction, products: Mol2dTuple
) -> MolTuple:
    """
    Recusrive sub-method for the method above.
    """
    single_step_product = generate_single_step_product(start_mol, reaction)
    if not single_step_product:
        return (start_mol,)
    else:
        for scenario in single_step_product:
            for p in scenario:
                products += _generate_multi_step_product_recursion(
                    p, reaction, products
                )

        return products


def _get_reactant_substructures(reaction: Reaction) -> list[Mol]:
    smarts: str = AllChem.ReactionToSmarts(reaction)
    substructures_str = smarts[: smarts.index(">")].split(".")
    substructures = [Chem.MolFromSmarts(s) for s in substructures_str]
    return substructures


def _subreaction_is_valid(start_mol: Mol, subreaction: Reaction) -> bool:
    substructures = _get_reactant_substructures(subreaction)
    for s in substructures:
        if start_mol.HasSubstructMatch(s):
            return True
    return False


def find_possible_reactions(
    start_mol: Mol,
    all_reactions: list[Reaction],
    solver_mode: bool,
) -> list[Reaction]:
    """
    Returns the reactions that can be performed on the given molecule.
    """
    possible_reactions = []
    for reaction in all_reactions:
        if solver_mode and reaction.num_reactants > 1:
            continue  # Solver mode will omit any reactions with more than one reactant

        for subreaction in reaction.reactions_list:
            if _subreaction_is_valid(start_mol, subreaction):
                possible_reactions.append(reaction)
                break

    return possible_reactions


def find_synthetic_pathway(
    start_mol: Mol,
    target_mol: Mol,
    all_reactions: list[Reaction],
    max_solver_steps: int,
) -> tuple[bool, list[Reaction], list[int]]:
    """
    Auto-solver, utilizes multi step products.
    Performs breadth first search (BFS) to find the SHORTEST possible reaction pathway to the target molecule.
    """
    if not target_mol:
        raise ValueError("Target molecule must be present to use the solver mode.")

    # Set of visited nodes to prevent loops
    visited = set()
    queue = deque()

    # Add the start_node to the queue and visited list
    queue.append(start_mol)
    visited.add(Chem.MolToSmiles(start_mol))

    # The format of this dictionary is:
    # key = product SMILES -> value = reactant SMILES
    parent: dict[str, str] = {}
    parent[Chem.MolToSmiles(start_mol)] = None  # startMol has no previous reactants

    # The format of this dictionary is:
    # key = (reactant SMILES, product SMILES) -> value = (reaction, product_index) that converts
    # the reactant into the product. Needed to reconstruct the reaction pathway later.
    choices: dict[tuple[str, str], tuple[Reaction, int]] = {}

    # BFS loop
    dist = 0
    path_found = False
    while queue and dist < max_solver_steps:
        for _ in range(len(queue)):
            cur = queue.popleft()
            if Chem.MolToSmiles(cur) == Chem.MolToSmiles(target_mol):
                path_found = True
                break

            possible_reactions = find_possible_reactions(
                cur, all_reactions, solver_mode=True
            )
            all_reactions_products: list[Mol2dTuple] = [
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

    # Reconstructing the reaction pathway
    rxn_pathway = []
    choice_pathway = []
    target_smiles = Chem.MolToSmiles(target_mol)

    if path_found:
        while parent[target_smiles] is not None:
            reaction, product_index = choices[(parent[target_smiles], target_smiles)]
            rxn_pathway.append(reaction)
            choice_pathway.append(product_index)
            target_smiles = parent[target_smiles]

        rxn_pathway.reverse()
        choice_pathway.reverse()

    return path_found, rxn_pathway, choice_pathway
