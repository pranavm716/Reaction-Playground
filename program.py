from rdkit.Chem import Draw
from rdkit.Chem.rdchem import Mol

from computations import (
    copy_mol,
    find_synthetic_pathway,
    generate_multi_step_product,
)
from PIL.Image import Image as PILImage

from datatypes import SolverModeImageData


# @dataclass
# class Program:
#     ui: UI
#     start_mol: Mol
#     target_mol: Mol | None
#     multi_step_react_mode: bool
#     max_num_solver_steps: int
#     frontend_images_folder: str
#     all_reactions: list[Reaction]
#     multiple_reactants_prompts: dict[str, list[str]]
#
#     def run_program(self):
#         proceed_to_playground_mode = True
#         if self.target_mol:
#             self.run_solver_mode()
#             # proceed_to_playground_mode = (
#             #     self.ui.get_user_input(
#             #         "\nWould you like to continue using Reaction Playground (y/n)? "
#             #     )
#             #     == "y"
#             # )
#             # if proceed_to_playground_mode:
#             #     self.ui.print("\nThis is your starting molecule:\n")
#
#         # if proceed_to_playground_mode:
#         #     self.run_playground_mode()
#
#     def run_solver_mode(self):
#         """
#         Runs the auto synthetic pathway solver.
#         """
#         image_file_paths: dict[str, list[str]] = {
#             "start_mol_src": [self.get_full_image_path("start_mol.png")],
#             "target_mol_src": [self.get_full_image_path("target_mol.png")],
#         }
#         self.ui.save_mol(self.start_mol, image_file_paths["start_mol_src"][0])
#         self.ui.save_mol(self.target_mol, image_file_paths["target_mol_src"][0])
#         path_found, reaction_pathway, choice_pathway = find_synthetic_pathway(
#             self.start_mol,
#             self.target_mol,
#             self.all_reactions,
#             self.max_num_solver_steps,
#         )
#
#         if not path_found:
#             return path_found, reaction_pathway, choice_pathway
#
#         current_mol = copy_mol(self.start_mol)
#         for step_number, (reaction, choice) in enumerate(
#             zip(reaction_pathway, choice_pathway), start=1
#         ):
#             products = generate_multi_step_product(current_mol, reaction)
#
#             step_full_image_path = self.get_full_image_path(f"{step_number - 1}.png")
#             self.ui.save_mol(current_mol, step_full_image_path)
#             image_file_paths[str(step_number - 1)] = [step_full_image_path]
#
#             num_products = len(products)
#             if num_products > 1:
#                 for i in range(1, num_products + 1):
#                     file_name = self.get_full_image_path(f"{step_number - 1}_{i}.png")
#                     image_file_paths[f"{step_number - 1}"].append(file_name)
#                     self.ui.save_mol(products[i - 1], file_name)
#
#             next_mol = products[choice]
#             current_mol = copy_mol(next_mol)
#
#         return path_found, reaction_pathway, choice_pathway, image_file_paths


# def get_full_image_path(file_name: str, frontend_images_folder: str) -> str:
#     return os.path.join(
#         os.path.dirname(__file__),
#         frontend_images_folder,
#         file_name,
#     )


def construct_mol_image(mol: Mol) -> PILImage:
    return Draw.MolToImage(mol)


# @dataclass
# class SolverModeImageData:
#     start_mol_image: PILImage
#     target_mol_image: PILImage
#     intermediate_steps_images


def run_solver_mode(
    start_mol: Mol, target_mol: Mol
) -> tuple[bool, int, list[str], list[int], PILImage, PILImage, SolverModeImageData]:
    """
    Runs the auto synthetic pathway solver.
    """

    start_mol_img = construct_mol_image(start_mol)
    target_mol_img = construct_mol_image(target_mol)
    solver_images: SolverModeImageData = []

    path_found, reaction_pathway, choice_pathway = find_synthetic_pathway(
        start_mol, target_mol
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

        product_image_info: list[tuple[PILImage, str]] = []
        for product_index, product in enumerate(products):
            product_img = construct_mol_image(product)
            product_img_filename = f"{step_number}_{product_index}.png"
            product_image_info.append((product_img, product_img_filename))

        solver_images.append((product_image_info, choice))
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
