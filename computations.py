from collections import deque

from config import Mol2dTuple
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol
from reaction import Reaction


def generate_unique_products(products: Mol2dTuple) -> Mol2dTuple:
    """
    Returns only the unique products in a 2d mol tuple. TODO: Make this less ugly
    """
    smiles = tuple(
        tuple(Chem.MolToSmiles(s) for s in scenarios) for scenarios in products
    )
    unique_smiles = tuple(set([tuple(set(x)) for x in smiles]))
    return tuple(
        tuple(Chem.MolFromSmiles(s) for s in scenarios) for scenarios in unique_smiles
    )


def generate_single_step_product(
    start_mol: Mol,
    reaction: Reaction,
    additional_reactants: tuple[Mol] | None = None,
):
    """
    Performs a reaction on a molecule once, and returns a 2d tuple of products.
    """
    all_reactants = (start_mol,)
    if additional_reactants:
        all_reactants += additional_reactants

    products: Mol2dTuple = tuple()
    for r in reaction.reactions_list:
        products += r.RunReactants(all_reactants)

    return generate_unique_products(products)


def generate_multi_step_product(start_mol: Mol, reaction: Reaction):
    """
    Recursively performs a reaction on a starting molecule and on its products until no more products can be formed.
    Returns the final tuple of products (always 1 x n).
    """
    products = _generate_multi_step_product_recursion(start_mol, reaction, tuple())
    return generate_unique_products((products,))[0]


def _generate_multi_step_product_recursion(
    start_mol: Mol, reaction: Reaction, products: Mol2dTuple
):
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


def find_possible_reactions(
    start_mol: Mol, all_rxns: list[Reaction], solver_mode: bool = False
):
    """
    Finds out which reactions are valid for the given molecule.
    """
    possible_reactions = []
    for obj in all_rxns:
        if solver_mode and obj.num_reactants > 1:
            continue  # Solver mode will omit any reactions with more than one reactant
        possible_scenarios = []
        for rxn in obj.reactions_list:
            rxn_smarts = AllChem.ReactionToSmarts(rxn)
            reactants_side = rxn_smarts[: rxn_smarts.index(">")]
            reactants_side = reactants_side.split(".")
            reaction_match = []
            for mol in reactants_side:
                reaction_match.append(
                    start_mol.HasSubstructMatch(Chem.MolFromSmarts(mol))
                )
            possible_scenarios.append(any(reaction_match))
        possible_reactions.append(any(possible_scenarios))
    possible_indices = [i for (i, val) in enumerate(possible_reactions) if val]
    return possible_indices


def find_synthetic_pathway(start_mol: Mol, target_mol: Mol, num_steps):
    """
    Auto-solver, utilizes multi step products.
    Performs breadth first search (BFS) to find the SHORTEST possible reaction pathway to the target molecule.
    """

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
    # key = (reactant, product) -> value = (rxn_index, choice_index) that converts the reactant into the product
    # Needed to reconstruct the reaction pathway later
    choices: dict[tuple[str, str], tuple[int, int]] = {}

    # BFS loop
    dist = 0
    path_found = False
    while queue and dist < num_steps:
        for _ in range(len(queue)):
            cur = queue.popleft()
            if Chem.MolToSmiles(cur) == Chem.MolToSmiles(target_mol):
                path_found = True
                break

            possible_indices = find_possible_reactions(cur, all_rxns, solver_mode=True)
            multi_step_products = [
                generate_multi_step_product(cur, rxn)
                for rxn in [all_rxns[i] for i in possible_indices]
            ]

            # Loop through all compatible reactions
            for i, msp in enumerate(multi_step_products):
                reaction_index = possible_indices[i]

                # Loop through all products of each compatible reaction
                for product_index, mol in enumerate(msp):
                    if Chem.MolToSmiles(mol) not in visited:
                        visited.add(Chem.MolToSmiles(mol))
                        queue.append(mol)

                        parent[Chem.MolToSmiles(mol)] = Chem.MolToSmiles(cur)
                        choices[
                            (Chem.MolToSmiles(cur), Chem.MolToSmiles(mol))
                        ] = (reaction_index, product_index)
        if path_found:
            break
        dist += 1

    # Reconstructing the reaction pathway
    rxn_pathway = []
    choice_pathway = []
    target_smiles = Chem.MolToSmiles(target_mol)
    if path_found:
        while parent[target_smiles] is not None:
            reaction_index, product_index = choices[
                (parent[target_smiles], target_smiles)
            ]
            rxn_pathway.append(reaction_index)
            choice_pathway.append(product_index)
            target_smiles = parent[target_smiles]

        rxn_pathway.reverse()
        choice_pathway.reverse()

    return path_found, rxn_pathway, choice_pathway
