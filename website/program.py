import copy

from PIL.Image import Image as PILImage
from rdkit.Chem.rdchem import Mol

from website.computations import (
    copy_mol,
    find_synthetic_pathway,
    generate_multi_step_product,
    get_reactant_position_of_mol_in_reaction,
)
from website.datatypes import SolverModeImageData, MolTuple
from website.fastapi_rdkit_utils import construct_mol_image
from website.reaction import Reaction


def run_solver_mode(
    start_mol: Mol, target_mol: Mol, all_reactions: list[Reaction]
) -> tuple[bool, int, list[str], list[int], PILImage, PILImage, SolverModeImageData]:
    """
    Runs the auto synthetic pathway solver. Returns information about the calculated synthetic pathway,
    including images of the molecules at each step.
    """

    start_mol_img = construct_mol_image(start_mol)
    target_mol_img = construct_mol_image(target_mol)

    # Includes the target molecule as the last entry, does not include the starting molecule
    solver_images: SolverModeImageData = {}

    path_found, reaction_pathway, choice_pathway = find_synthetic_pathway(
        start_mol, target_mol, all_reactions
    )
    num_steps = len(choice_pathway)

    # Reconstruct the synthetic pathway
    reaction_names: list[str] = []
    current_mol = copy_mol(start_mol)
    for step_number, (reaction, choice) in enumerate(
        zip(reaction_pathway, choice_pathway)
    ):
        reaction_names.append(reaction.name)
        products = generate_multi_step_product(current_mol, reaction)

        product_images = [construct_mol_image(product) for product in products]
        solver_images[step_number] = product_images

        current_mol = products[choice]

    return (
        path_found,
        num_steps,
        reaction_names,
        choice_pathway,
        start_mol_img,
        target_mol_img,
        solver_images,
    )


def get_missing_reactant_prompts(current_mol: Mol, reaction: Reaction) -> list[str]:
    """
    This method will return the prompts for reactions that require additional reactants.
    """

    reactant_position = get_reactant_position_of_mol_in_reaction(current_mol, reaction)
    multiple_reactant_prompts = copy.copy(reaction.multiple_reactants_prompts)
    multiple_reactant_prompts.pop(reactant_position)
    return multiple_reactant_prompts
