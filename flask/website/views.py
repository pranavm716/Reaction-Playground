from flask import Blueprint, flash, redirect, render_template, request, session, url_for
from rdkit import Chem

from config import read_config
from program import Program
from ui import CLI

views = Blueprint("views", __name__)


@views.route("/")
def home():
    return render_template("home.jinja")


@views.route("/", methods=["POST"])
def process_init():
    start_mol_SMILES = request.form.get("start_mol_SMILES")
    target_mol_SMILES = request.form.get("target_mol_SMILES")

    print(start_mol_SMILES, target_mol_SMILES)
    try:
        start_mol = Chem.MolFromSmiles(start_mol_SMILES)
        Chem.SanitizeMol(start_mol)
        target_mol = Chem.MolFromSmiles(target_mol_SMILES)
        Chem.SanitizeMol(target_mol)
    except Exception as e:
        flash("Invalid SMILES. Please enter valid SMILES.", category="error")
        return render_template("home.jinja")

    config = read_config("./config.json")

    if config.disable_rdkit_warnings:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.warning")

    program = Program(
        ui=CLI(),
        start_mol=start_mol,
        target_mol=target_mol,
        multi_step_react_mode=config.multi_step_react_mode,
        max_num_solver_steps=config.max_num_solver_steps,
        all_reactions_file_path=config.all_reactions_file_path,
        multiple_reactants_prompts=config.multiple_reactants_prompts,
    )
    session["program"] = program
    return redirect(url_for("views.program"))


@views.route("/program")
def program():
    program = session["program"]
    return render_template(
        "program.jinja",
        start_mol_SMILES=Chem.MolToSMILES(program.start_mol),
        target_mol_SMILES=Chem.MolToSMILES(program.target_mol),
    )
