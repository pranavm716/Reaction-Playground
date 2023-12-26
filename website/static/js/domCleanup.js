function cleanUpDisplayReactions() {
    document.getElementById("reactionForm").remove();
    document.getElementById("reactionSummary").remove();
}

function cleanUpMultipleProducts(productIndex) {
    highlightSelectedProductAndDisableOnClicks(productIndex);
    removeMultipleProductsPrompt();
}

function highlightSelectedProductAndDisableOnClicks(productIndex) {
    let productsDiv = document.getElementById("productsDiv");

    // Select images inside the div with id 'productsDiv' whose ids start with 'product_'
    const products = $('#productsDiv img[id^="product_"]');
    for (let index = 0; index < products.length; index++) {
        let product = products[index];
        product.removeAttribute("onclick");
        if (index === productIndex) {
            product.className = "card-img-top bordered-image"
        }
    }
    productsDiv.removeAttribute("id");
}

function removeMultipleProductsPrompt() {
    document.getElementById("multipleProductsPrompt").remove();
}

function updateCurrentMolDiv() {
    document.getElementById("currentMolDiv").removeAttribute("id");
}