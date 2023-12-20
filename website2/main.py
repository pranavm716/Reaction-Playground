import pathlib
from typing import Annotated

from fastapi import FastAPI, Request, Form, HTTPException
from starlette.responses import HTMLResponse, RedirectResponse
from starlette.staticfiles import StaticFiles
from starlette.templating import Jinja2Templates

from program import run_solver_mode
from reaction import read_all_reactions_from_file
from website2.fastapi_rdkit_utils import (
    start_and_target_mols_are_valid,
    construct_query_url,
    img_to_base64,
    get_mol_from_smiles,
)
import tempfile
import os


MULTI_STEP_REACT_MODE = True
MAX_NUM_SOLVER_STEPS = 15
DISABLE_RDKIT_WARNINGS = True

if DISABLE_RDKIT_WARNINGS:
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.warning")

ALL_REACTIONS_FILE_PATH = pathlib.Path("./all_reactions.yaml")
ALL_REACTIONS = read_all_reactions_from_file(ALL_REACTIONS_FILE_PATH)

app = FastAPI()

static_dir = pathlib.Path(__file__).parent.resolve() / "static"
templates_dir = pathlib.Path(__file__).parent.resolve() / "templates"

app.mount("/static", StaticFiles(directory=static_dir), name="static")
templates = Jinja2Templates(directory=templates_dir)


@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    return templates.TemplateResponse("home.jinja", {"request": request})


@app.post("/", response_class=HTMLResponse)
async def run_program(
    request: Request,
    start_mol_smiles: Annotated[str, Form()],
    target_mol_smiles: Annotated[str, Form()],
):
    if not start_and_target_mols_are_valid(start_mol_smiles, target_mol_smiles):
        raise HTTPException(status_code=400, detail="Invalid form data")

    if target_mol_smiles:
        solver_mode_url = construct_query_url(
            app,
            "solver_mode",
            start_mol_smiles=start_mol_smiles,
            target_mol_smiles=target_mol_smiles,
        )
        return RedirectResponse(solver_mode_url, status_code=303)


@app.get("/solver-mode", response_class=HTMLResponse)
async def solver_mode(request: Request, start_mol_smiles: str, target_mol_smiles: str):
    start_mol = get_mol_from_smiles(start_mol_smiles)
    target_mol = get_mol_from_smiles(target_mol_smiles)
    (
        path_found,
        num_steps,
        reaction_names,
        choice_pathway,
        start_mol_img,
        target_mol_img,
        solver_images,
    ) = run_solver_mode(start_mol, target_mol)

    with tempfile.TemporaryDirectory() as solver_images_dir:
        pass
        # start_mol_img.save(pathlib.Path(solver_images_dir) / "start_mol.png", "PNG")
        # target_mol_img.save(pathlib.Path(solver_images_dir) / "target_mol.png", "PNG")
        # for imgs, _ in solver_images:
        #     for img, fp in imgs:
        #         img.save(pathlib.Path(solver_images_dir) / fp, "PNG")

        # print(os.listdir(solver_images_dir))
        # return templates.TemplateResponse(
        #     "solver_mode.jinja",
        #     {
        #         "request": request,
        #         "path_found": path_found,
        #         "num_steps": num_steps,
        #         "reaction_names": reaction_names,
        #         "solver_images_dir": pathlib.Path(solver_images_dir).resolve(),
        #         "max_num_solver_steps": MAX_NUM_SOLVER_STEPS,
        #     },
        # )
    return templates.TemplateResponse(
        "solver_mode.jinja",
        {
            "request": request,
            "path_found": path_found,
            "num_steps": num_steps,
            "reaction_names": reaction_names,
            "max_num_solver_steps": MAX_NUM_SOLVER_STEPS,
            "start_mol": img_to_base64(start_mol_img),
        },
    )
