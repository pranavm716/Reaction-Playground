from typing import Iterable
from typing import Annotated

from fastapi import FastAPI, Request, Form, HTTPException
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from starlette.responses import HTMLResponse, RedirectResponse
from starlette.staticfiles import StaticFiles
from starlette.templating import Jinja2Templates

from program import run_solver_mode, get_missing_reactant_prompts
from website.computations import (
    find_possible_reaction_keys,
    generate_single_step_product,
    generate_multi_step_product,
    copy_mol,
    get_reactant_position_of_mol_in_reaction,
    ALL_REACTIONS,
    get_reactions_from_keys,
)
from website.config import (
    MAX_NUM_SOLVER_STEPS,
    MULTI_STEP_REACT_MODE,
    STATIC_DIR,
    TEMPLATES_DIR,
)
from website.datatypes import Mol2dTuple
from website.fastapi_rdkit_utils import (
    start_and_target_mols_are_valid,
    construct_query_url,
    img_to_base64,
    mol_to_base64,
)
from website.reaction import Reaction, ReactionKey

app = FastAPI()


# Mount static folder
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

# Create jinja2 template engine
templates = Jinja2Templates(directory=TEMPLATES_DIR)


@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    """Allows the user to enter SMILES for a starting and a target molecule."""
    return templates.TemplateResponse("home.jinja", {"request": request})


@app.post("/", response_class=RedirectResponse)
async def run_program(
    request: Request,
    start_mol_smiles: Annotated[str, Form()],
    target_mol_smiles: Annotated[str | None, Form()] = None,
):
    """
    Runs the solver mode if start and target SMILES are present and playground mode if only start SMILES is present.
    Also validates the SMILES and raises an HTTPException if invalid.
    """

    if not start_and_target_mols_are_valid(start_mol_smiles, target_mol_smiles):
        raise HTTPException(
            status_code=400, detail="Invalid starting and/or target molecule SMILES."
        )

    if target_mol_smiles:
        solver_mode_url = construct_query_url(
            app,
            "solver_mode",
            start_mol_smiles=start_mol_smiles,
            target_mol_smiles=target_mol_smiles,
        )
        return RedirectResponse(solver_mode_url, status_code=303)
    else:
        playground_mode_url = construct_query_url(
            app, "playground_mode", mol_smiles=start_mol_smiles
        )
        return RedirectResponse(playground_mode_url, status_code=303)


@app.get("/solver-mode", response_class=HTMLResponse)
async def solver_mode(request: Request, start_mol_smiles: str, target_mol_smiles: str):
    """Displays the steps needed to synthesize the target molecule from the start molecule."""

    start_mol = Chem.MolFromSmiles(start_mol_smiles)
    target_mol = Chem.MolFromSmiles(target_mol_smiles)

    (
        path_found,
        num_steps,
        reaction_names,
        choice_pathway,
        start_mol_img,
        target_mol_img,
        solver_images,
    ) = run_solver_mode(start_mol, target_mol)

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


@app.get("/playground-mode", response_class=HTMLResponse)
async def playground_mode(request: Request, mol_smiles: str):
    """
    Starts the playground mode loop, where users can experiment with applying reactions freely to molecules.
    """

    current_mol: Mol = Chem.MolFromSmiles(mol_smiles)
    app.state.current_mol = current_mol

    return templates.TemplateResponse("playground_mode.jinja", {"request": request})


@app.get("/playground-mode/choose-reaction", response_class=HTMLResponse)
async def playground_mode_choose_reaction(request: Request):
    """Displays the valid reactions for a molecule."""

    current_mol: Mol = app.state.current_mol

    possible_reaction_keys = find_possible_reaction_keys(current_mol, solver_mode=False)
    app.state.possible_reaction_keys = possible_reaction_keys

    return templates.TemplateResponse(
        "playground_mode_choose_reaction.jinja",
        {
            "request": request,
            "current_mol": mol_to_base64(current_mol),
            "current_mol_smiles": Chem.MolToSmiles(current_mol),
            "possible_reactions": get_reactions_from_keys(possible_reaction_keys),
        },
    )


@app.post("/playground-mode/display-products", response_class=HTMLResponse)
async def playground_mode_display_products(
    request: Request, choice: Annotated[str, Form()]
):
    """Displays the products that result from running the chosen reaction."""

    current_mol: Mol = app.state.current_mol
    possible_reaction_keys: list[ReactionKey] = app.state.possible_reaction_keys

    if choice in {str(i) for i in range(len(possible_reaction_keys))}:
        chosen_reaction_key = possible_reaction_keys[int(choice)]
        app.state.chosen_reaction_key = chosen_reaction_key
    else:
        raise HTTPException(
            status_code=404, detail=f"Invalid reaction choice {choice!r}."
        )

    chosen_reaction = ALL_REACTIONS[chosen_reaction_key]
    if chosen_reaction.num_reactants > 1:
        # Handling the reactions that require additional reactants (need to prompt user)
        prompts = get_missing_reactant_prompts(current_mol, chosen_reaction_key)
        return templates.TemplateResponse(
            "playground_mode_add_reactants.jinja",
            {"request": request, "prompts": prompts},
        )
    elif not MULTI_STEP_REACT_MODE:
        products = generate_single_step_product(current_mol, chosen_reaction_key)
    else:
        products = (generate_multi_step_product(current_mol, chosen_reaction_key),)
    app.state.products = products

    return templates.TemplateResponse(
        "playground_mode_display_products.jinja",
        {
            "request": request,
            "chosen_reaction": chosen_reaction,
            "products": [
                [mol_to_base64(product) for product in scenario]
                for scenario in products
            ],
            "product_smiles": [
                [Chem.MolToSmiles(product) for product in scenario]
                for scenario in products
            ],
        },
    )


@app.post("/playground-mode/process-added-reactants", response_class=HTMLResponse)
def playground_mode_process_added_reactants(
    request: Request, extra_reactant_smiles: Annotated[list[str], Form()]
):
    """
    Takes the SMILES for the extra needed reactants (for reactions with multiple products) and runs
    the chosen reaction.
    """

    # TODO: single step products are still not supported in the UI
    current_mol: Mol = app.state.current_mol
    chosen_reaction_key: ReactionKey = app.state.chosen_reaction_key

    reactants: list[Mol] = [
        Chem.MolFromSmiles(smiles) for smiles in extra_reactant_smiles
    ]
    reactant_position = get_reactant_position_of_mol_in_reaction(
        current_mol, chosen_reaction_key
    )
    reactants.insert(reactant_position, current_mol)

    products = generate_single_step_product(tuple(reactants), chosen_reaction_key)
    app.state.products = products

    return templates.TemplateResponse(
        "playground_mode_display_products.jinja",
        {
            "request": request,
            "chosen_reaction": ALL_REACTIONS[chosen_reaction_key],
            "products": [
                [mol_to_base64(product) for product in scenario]
                for scenario in products
            ],
            "product_smiles": [
                [Chem.MolToSmiles(product) for product in scenario]
                for scenario in products
            ],
        },
    )


@app.post("/playground-mode/choose-product", response_class=HTMLResponse)
def playground_mode_choose_product(
    request: Request, product_index: Annotated[str, Form()]
):
    """Updates the app's state with the chosen product."""

    # Assume multi step react mode for now
    products: Mol2dTuple = app.state.products
    if product_index not in {str(i) for i in range(len(products[0]))}:
        raise HTTPException(
            status_code=404, detail=f"Invalid product index {product_index!r}"
        )

    product_index = int(product_index)
    scenario_index = 0

    next_mol = products[scenario_index][product_index]
    app.state.current_mol = copy_mol(next_mol)


@app.get("/reactions/all", response_model=Iterable[Reaction])
def all_reactions(request: Request):
    return ALL_REACTIONS.values()
