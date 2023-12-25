// Displays the valid reactions for a molecule.
// Called from playground_mode.jinja
async function displayReactions() {
    const response = await fetch("/playground-mode/choose-reaction", {
        method: "GET",
    });
    updateContent(await response.text());
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
        await chooseProduct(0);
    }
}

// Updates the app's state with the chosen product.
// Called from reaction_playground_display_products.jinja
async function chooseProduct(productIndex) {
    disableOnClickForAllProducts();
    await fetch("/playground-mode/choose-product", {
        method: "POST",
        body: new URLSearchParams({product_index: productIndex}),
    });
    await displayReactions();
}
