// Displays the valid reactions for a molecule.
// Called from playground_mode.jinja
async function displayReactions(onlyOneProductPreviouslyDisplayed = true) {
    const response = await fetch("/playground-mode/choose-reaction", {
        method: "GET",
    });
    updateContent(await response.text());

    if (!onlyOneProductPreviouslyDisplayed) {
        document.getElementById("currentMolDiv").remove();
    } else {
        updateCurrentMolDiv();
    }
}

// Displays the products that result from running the chosen reaction
// Called from reaction_playground_choose_reaction.jinja
async function displayProducts(choice) {
    cleanUpDisplayReactions()
    const response = await fetch("/playground-mode/display-products", {
        method: "POST",
        body: new URLSearchParams({choice: choice}),
    });
    updateContent(await response.text());

    if (onlyOneProductDisplayed()) {
        document.getElementById("productsDiv").remove();
        await chooseProduct(0, true);
    }
}

// Updates the app's state with the chosen product.
// Called from reaction_playground_display_products.jinja
async function chooseProduct(productIndex, onlyOneProductDisplayed = false) {
    if (!onlyOneProductDisplayed) {
        cleanUpMultipleProducts(productIndex);
    }

    await fetch("/playground-mode/choose-product", {
        method: "POST",
        body: new URLSearchParams({product_index: productIndex}),
    });
    await displayReactions(onlyOneProductDisplayed);
}

async function processAddedReactants(numSmiles) {
    let smilesList = [];

    for (let index = 0; index < numSmiles; index++) {
        const inputId = "extra_reactant_smiles_" + index;
        const inputValue = document.getElementById(inputId).value;
        smilesList.push(inputValue);
    }

    const response = await fetch("/playground-mode/process-added-reactants", {
        method: "POST",
        body: new URLSearchParams({extra_reactant_smiles: smilesList}),
    });
    updateContent(await response.text());

    cleanUpExtraReactantsDiv();
    if (onlyOneProductDisplayed()) {
        document.getElementById("productsDiv").remove();
        await chooseProduct(0, true);
    }
}
