from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol

from computations import find_synthetic_pathway
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

    def run_program(self):
        solver_mode = bool(self.target_mol)
        proceed_to_playground_mode = True
        if solver_mode:
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
        history = [self.start_mol]
        while True:
            pass
        # TODO: rest of this method


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
