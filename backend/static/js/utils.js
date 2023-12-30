function updateContent(html) {
    document.getElementById("parentContainer").innerHTML += html;
}

function onlyOneProductDisplayed() {
    let productsDiv = document.getElementById("productsDiv");
    if (!productsDiv) {
        return false;
    }

    // Select images inside the div with id 'productsDiv' whose ids start with 'product_'
    const numberOfImages = $('#productsDiv img[id^="product_"]').length;
    return numberOfImages === 1;
}
