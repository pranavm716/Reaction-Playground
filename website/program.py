from PIL.Image import Image as PILImage
from rdkit.Chem.rdchem import Mol

from website.computations import (
    copy_mol,
    find_synthetic_pathway,
    generate_multi_step_product,
)
from website.datatypes import SolverModeImageData
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

    # def run_playground_mode(self) -> None:
    #     """
    #     Runs the playground mode loop, where users can experiment with applying reactions
    #     freely to molecules.
    #     """
    #     history: list[Mol] = []
    #     current_mol = copy_mol(self.start_mol)
    #     while True:
    #         possible_reactions = find_possible_reactions(
    #             current_mol, self.all_reactions, solver_mode=False
    #         )
    #         self.ui.display_compatible_reactions(current_mol, possible_reactions)

    #         if history:
    #             self.ui.print("[b] Back\n")
    #         self.ui.print("[q] Quit\n")

    #         choice = self.ui.get_user_input()
    #         if choice == "q":
    #             break
    #         elif history and choice == "b":
    #             current_mol = history.pop()
    #             continue
    #         elif choice in [str(i) for i in range(1, len(possible_reactions) + 1)]:
    #             chosen_reaction = possible_reactions[int(choice) - 1]
    #         else:
    #             self.ui.print("Invalid input!")
    #             continue

    #         if chosen_reaction.num_reactants > 1:
    #             # Handling the reactions that require additional reactants (need to prompt user)
    #             reactants = self.get_missing_reactants(current_mol, chosen_reaction)
    #             products = generate_single_step_product(reactants, chosen_reaction)
    #         elif not self.multi_step_react_mode:
    #             products = generate_single_step_product(current_mol, chosen_reaction)
    #         else:
    #             products = (generate_multi_step_product(current_mol, chosen_reaction),)

    #         if not products:
    #             self.ui.print("\nNo reaction!\n")
    #             continue

    #         self.ui.print()
    #         if len(products[0]) == 1 and len(products) == 1:
    #             self.ui.print("The product is:")
    #         else:
    #             self.ui.print("The products are:")

    #         chosen_scenario, chosen_product = 0, 0
    #         if len(products) > 1:
    #             self.ui.display_2d_mol_tuple(products)
    #             chosen_scenario = self.ui.prompt_for_scenario_index(products)
    #         if len(products[0]) > 1:
    #             self.ui.save_mol_tuple(products[chosen_scenario])
    #             chosen_product = self.ui.prompt_for_product_index(
    #                 products[chosen_scenario]
    #             )

    #         history.append(current_mol)
    #         next_mol = products[chosen_scenario][chosen_product]
    #         current_mol = copy_mol(next_mol)

    # def get_missing_reactants(self, current_mol: Mol, reaction: Reaction) -> MolTuple:
    #     """
    #     For reactions that require additional reactants, this method will prompt the
    #     user for those additional missing reactants.
    #     """
    #     prompts = self.multiple_reactants_prompts[reaction.name]
    #     reactant_position = get_reactant_position_of_mol_in_reaction(
    #         current_mol, reaction
    #     )

    #     reactants: list[Mol] = []
    #     for i, prompt in enumerate(prompts):
    #         if i == reactant_position:
    #             reactants.append(current_mol)
    #         else:
    #             missing_reactant_smiles = self.ui.get_user_input(prompt)
    #             missing_reactant = Chem.MolFromSmiles(missing_reactant_smiles)

    #             self.ui.print("\nYou entered:")
    #             self.ui.save_mol(missing_reactant)
    #             self.ui.print()

    #             reactants.append(missing_reactant)

    #     return tuple(reactants)

    # def get_full_image_path(self, file_name: str) -> str:
    #     return os.path.join(
    #         os.path.dirname(__file__),
    #         self.frontend_images_folder,
    #         file_name,
    #     )
