import os

from flask import Blueprint, flash, redirect, render_template, request, session, url_for
from rdkit import Chem

from backend.config import read_config
from backend.program import Program
from backend.ui import CLI

views = Blueprint("views", __name__)


@views.route("/", methods=["GET", "POST"])
def home():
    def flash_invalid_mol(mol_type: str):
        flash(
            f"Invalid SMILES for {mol_type} molecule. Please enter valid SMILES.",
            category="error",
        )

    if request.method == "GET":
        print("GET")
        return render_template("home.jinja")

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

    start_mol_smiles = session.get("start_mol_smiles")
    start_mol = Chem.MolFromSmiles(start_mol_smiles)

    target_mol_smiles = session.get("target_mol_smiles")
    target_mol = Chem.MolFromSmiles(target_mol_smiles) if target_mol_smiles else None
    
    program = Program(
        ui=CLI(),
        start_mol=start_mol,
        target_mol=target_mol,
        multi_step_react_mode=config.multi_step_react_mode,
        max_num_solver_steps=config.max_num_solver_steps,
        all_reactions_file_path=config.all_reactions_file_path,
        multiple_reactants_prompts=config.multiple_reactants_prompts,
    )

    return render_template(
        "solver_mode.jinja",
        start_mol_SMILES=start_mol_smiles,
        target_mol_SMILES=target_mol_smiles,
    )
