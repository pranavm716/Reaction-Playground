import os

from flask import Blueprint, flash, redirect, render_template, request, session, url_for
from rdkit import Chem

from backend.config import read_config
from backend.interfaces.flask_ui import FlaskUI
from backend.program import Program
from backend.reaction import Reaction

views = Blueprint("views", __name__)
DYNAMIC_IMAGES_DIR = os.path.join(os.path.dirname(__file__), "./static/images/dynamic")


def clear_directory(directory):
    for filename in os.listdir(directory):
        if filename == ".gitignore":
            continue
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)


@views.route("/", methods=["GET"])
def home():
    return render_template("home.jinja")


@views.route("/", methods=["POST"])
def store_start_and_target_molecules():
    def flash_invalid_mol(mol_type: str):
        flash(
            f"Invalid SMILES for {mol_type} molecule. Please enter valid SMILES.",
            category="error",
        )

    clear_directory(DYNAMIC_IMAGES_DIR)
    start_mol_smiles = request.form.get("start_mol_smiles")
    target_mol_smiles = request.form.get("target_mol_smiles")
    print(f"{start_mol_smiles=}", f"{target_mol_smiles=}")

    invalid_string = ""
    if not start_mol_smiles:
        invalid_string += "starting"

    try:
        start_mol = Chem.MolFromSmiles(start_mol_smiles)
        Chem.SanitizeMol(start_mol)
    except Exception:
        invalid_string += "starting"

    try:
        target_mol = Chem.MolFromSmiles(target_mol_smiles)
        Chem.SanitizeMol(target_mol)
    except Exception:
        invalid_string += " and target"

    if invalid_string:
        flash_invalid_mol(invalid_string)
        return redirect(request.referrer)
    else:
        session["start_mol_smiles"] = start_mol_smiles
        if target_mol_smiles:
            session["target_mol_smiles"] = target_mol_smiles
        return redirect(url_for("views.run_program"))


@views.route("/run-program", methods=["GET"])
def run_program():
    config_file_location = os.path.join(
        os.path.dirname(__file__), "../../backend/config.json"
    )
    config = read_config(config_file_location)
    config.all_reactions_file_path = os.path.join(
        os.path.dirname(__file__), "../../backend", config.all_reactions_file_path
    )

    if config.disable_rdkit_warnings:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.warning")

    session["max_num_solver_steps"] = config.max_num_solver_steps

    start_mol_smiles = session["start_mol_smiles"]
    start_mol = Chem.MolFromSmiles(start_mol_smiles)

    target_mol_smiles = session["target_mol_smiles"]
    target_mol = Chem.MolFromSmiles(target_mol_smiles) if target_mol_smiles else None

    program = Program(
        ui=FlaskUI(),
        start_mol=start_mol,
        target_mol=target_mol,
        multi_step_react_mode=config.multi_step_react_mode,
        max_num_solver_steps=config.max_num_solver_steps,
        frontend_images_folder=config.frontend_images_folder,
        all_reactions_file_path=config.all_reactions_file_path,
        multiple_reactants_prompts=config.multiple_reactants_prompts,
    )

    if target_mol:
        (
            path_found,
            reaction_pathway,
            choice_pathway,
            image_file_paths,
        ) = program.run_solver_mode()
        return run_solver_mode(
            path_found, reaction_pathway, choice_pathway, image_file_paths
        )
    else:
        pass
        # program.run_playground_mode()
        # return run_playground_mode(program)


def run_solver_mode(
    path_found: bool,
    reaction_pathway: list[Reaction],
    choice_pathway: list[int],
    image_file_paths: dict[str, list[str]],
):
    
    session["path_found"] = path_found
    session["choice_pathway"] = choice_pathway
    session["reaction_names"] = [reaction.name for reaction in reaction_pathway]

    relative_image_file_paths = {}
    for k, v in image_file_paths.items():
        relative_image_file_paths[k] = [
            os.path.relpath(file_path, os.path.dirname(__file__)) for file_path in v
        ]
    session.update(relative_image_file_paths)
    print(f"{relative_image_file_paths=}")

    return render_template("solver_mode.jinja")


# def run_playground_mode(program: Program) -> Response:
#     pass
