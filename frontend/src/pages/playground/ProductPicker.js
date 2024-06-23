import React, { useEffect } from "react";
import MolImage from "../../components/MolImage";
import {ArrowWithReactionName, plusIcon} from "../../components/SmallUIComponents";

const ProductPicker = ({ products, setSmiles, reactionName, molImage, missingReactantSmilesPicked, missingReactantEncodings }) => {
    let reactantMolRow;
    if (missingReactantSmilesPicked) {
        reactantMolRow = <div className="mol-row">
            {molImage}
            {missingReactantSmilesPicked.map((smiles, index) => (
                <React.Fragment key={smiles}>
                    {plusIcon}
                    <MolImage smiles={smiles} encoding={missingReactantEncodings[index]} />
                </React.Fragment>
            ))}
        </div>
    } else {
        reactantMolRow = molImage;
    }

    useEffect(() => {
        // If there is only one product, there is no choice the user has to make
        if (products.length === 1) {
            setSmiles(products[0].smiles);
        }
    }, [products, setSmiles]);

    return (
        <>
            {reactantMolRow}
            <ArrowWithReactionName reactionName={reactionName} />

            {/* list of clickable products */}
            <div className="mol-row">
                {products.length > 1 && products.map((product, index) => (
                    <React.Fragment key={product.smiles}>
                        <MolImage smiles={product.smiles} encoding={product.encoding} setSmiles={setSmiles} />
                        {index < products.length - 1 && plusIcon}
                    </React.Fragment>
                ))}
            </div>
            <p className="multi-product-text">
                This produces {products.length} products. Click on the one you want to analyze next.
            </p>
        </>
    );
}

export default ProductPicker;