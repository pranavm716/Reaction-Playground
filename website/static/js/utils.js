function updateContent(html) {
    document.getElementById("parentContainer").innerHTML += html
}

function onlyOneProductDisplayed() {
    let result = false;
    if ($('#productsDiv').length) {
        // Select images inside the div with id 'productsDiv' whose ids start with 'product_'
        const numberOfImages = $('#productsDiv img[id^="product_"]').length;
        result = numberOfImages === 1;
    }
    return result;
}