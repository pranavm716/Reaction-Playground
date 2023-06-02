from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol

from computations import (
    copy_mol,
    find_possible_reactions,
    find_synthetic_pathway,
    generate_multi_step_product,
    generate_single_step_product,
)
from reaction import Reaction
from ui import UI


class Program:
    def __init__(
        self,
        ui: UI,
        start_mol: Mol,
        target_mol: Mol | None,
        multi_step_react_mode: bool,
        max_num_solver_steps: int,
        all_reactions_file_path: str,
    ) -> None:
        self.ui = ui
        self.start_mol = start_mol
        self.target_mol = target_mol
        self.multi_step_react_mode = multi_step_react_mode
        self.max_num_solver_steps = max_num_solver_steps
        self.all_reactions = read_all_reactions_from_file(all_reactions_file_path)
        self.solver_mode = bool(target_mol)

    def run_program(self):
        proceed_to_playground_mode = True
        if self.solver_mode:
            proceed_to_playground_mode = self._run_solver_mode()

        if not proceed_to_playground_mode:
            return

        self._run_playground_mode()

    def _run_solver_mode(self) -> bool:
        """
        Runs the auto synthetic pathway solver.
        Returns whether or not the user would like to continue using the tool in playground mode.
        """
        path_found, reaction_pathway, choice_pathway = find_synthetic_pathway(
            self.start_mol,
            self.target_mol,
            self.all_reactions,
            self.max_num_solver_steps,
        )
        self.ui.display_solver_mode_intro(self.start_mol, self.target_mol)
        # TODO: rest of this method
        return False

    def _run_playground_mode(self):
        history = []
        current_mol = copy_mol(self.start_mol)
        while True:
            possible_reactions = find_possible_reactions(
                current_mol, self.all_reactions, self.solver_mode
            )
            self.ui.display_compatible_reactions(current_mol, possible_reactions)

            if history:
                self.ui.print("[b] Back\n")
            self.ui.print("[q] Quit\n")

            choice = self.ui.get_user_input()
            if choice == "q":
                break
            elif history and choice == "b":
                current_mol = history.pop()
                continue
            elif choice in [str(i) for i in range(1, len(possible_reactions) + 1)]:
                chosen_reaction = possible_reactions[int(choice) - 1]
            else:
                self.ui.print("Invalid input!")
                continue

            if chosen_reaction.num_reactants > 1:
                # TODO: Handling the reactions that require additional reactants (need to prompt user)
                pass
            elif not self.multi_step_react_mode:
                products = generate_single_step_product(current_mol, chosen_reaction)
            else:
                products = (generate_multi_step_product(current_mol, chosen_reaction),)

            if not products:
                self.ui.print("\nNo reaction!\n")
                continue

            self.ui.print()
            if len(products[0]) == 1 and len(products) == 1:
                self.ui.print("The product is:")
            else:
                self.ui.print("The products are:")

            chosen_scenario, chosen_product = 0, 0
            if len(products) > 1:
                self.ui.display_2d_mol_tuple(products)
                chosen_scenario = self.ui.prompt_for_scenario_index(products)
            if len(products[0]) > 1:
                self.ui.display_mol_tuple(products[chosen_scenario])
                chosen_product = self.ui.prompt_for_product_index(
                    products[chosen_scenario]
                )

            current_mol = copy_mol(products[chosen_scenario][chosen_product])
            history.append(current_mol)


def read_all_reactions_from_file(all_reactions_file_path: str) -> list[Reaction]:
    """
    Returns a list of all available reactions.
    """

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
