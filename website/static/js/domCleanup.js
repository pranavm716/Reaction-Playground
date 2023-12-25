function cleanUpDisplayReactions() {
    document.getElementById("reactionForm").remove();
    document.getElementById("reactionSummary").remove();
}

function disableOnClickForAllProducts() {
    let productsDiv = document.getElementById("productsDiv");
    // Select images inside the div with id 'productsDiv' whose ids start with 'product_'
    const products = $('#productsDiv img[id^="product_"]');
    for (let product of products) {
        product.removeAttribute("onclick");
    }
    productsDiv.removeAttribute("id");
}
