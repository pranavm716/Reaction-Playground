import pathlib
from contextlib import asynccontextmanager
from typing import Annotated

from fastapi import FastAPI, Request, Form, HTTPException
from starlette.responses import HTMLResponse, RedirectResponse
from starlette.staticfiles import StaticFiles
from starlette.templating import Jinja2Templates

from program import run_solver_mode
from website.config import (
    MAX_NUM_SOLVER_STEPS,
    DISABLE_RDKIT_WARNINGS,
    ALL_REACTIONS_FILE_PATH,
)
from website.fastapi_rdkit_utils import (
    start_and_target_mols_are_valid,
    construct_query_url,
    img_to_base64,
    get_mol_from_smiles,
)
from website.reaction import read_all_reactions_from_file


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Disables rdkit warnings if configured and loads all reactions from a file."""

    if DISABLE_RDKIT_WARNINGS:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.warning")

    all_reactions = read_all_reactions_from_file(ALL_REACTIONS_FILE_PATH)
    app.state.all_reactions = all_reactions  # noqa

    yield


app = FastAPI(lifespan=lifespan)


# Mount static folder
static_dir = pathlib.Path(__file__).parent.resolve() / "static"
app.mount("/static", StaticFiles(directory=static_dir), name="static")

# Create jinja2 template engine
templates_dir = pathlib.Path(__file__).parent.resolve() / "templates"
templates = Jinja2Templates(directory=templates_dir)


@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    """Allows the user to enter SMILES for a starting and a target molecule."""
    return templates.TemplateResponse("home.jinja", {"request": request})


@app.post("/", response_class=RedirectResponse)
async def run_program(
    request: Request,
    start_mol_smiles: Annotated[str, Form()],
    target_mol_smiles: Annotated[str, Form()],
):
    """
    Runs the solver mode if start and target SMILES are present and playground mode if only start SMILES is present.
    Also validates the SMILES and raises an HTTPException if invalid.
    """
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
    """Displays the steps needed to synthesize the target molecule from the start molecule."""

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
    ) = run_solver_mode(
        start_mol=start_mol,
        target_mol=target_mol,
        all_reactions=app.state.all_reactions,  # noqa
    )

    solver_run_metrics = {
        "path_found": path_found,
        "num_steps": num_steps,
        "reaction_names": reaction_names,
        "choice_pathway": choice_pathway,
    }

    image_encodings = {
        "start_mol": img_to_base64(start_mol_img),
        "target_mol": img_to_base64(target_mol_img),
    }
    for step_number, product_images in solver_images.items():
        image_encodings[str(step_number)] = [
            img_to_base64(img) for img in product_images
        ]

    return templates.TemplateResponse(
        "solver_mode.jinja",
        {
            "request": request,
            **solver_run_metrics,
            "image_encodings": image_encodings,
            "max_num_solver_steps": MAX_NUM_SOLVER_STEPS,
        },
    )
