import React, { useEffect } from "react";
import MolImage from "./MolImage";
import {ArrowWithReactionName, plusIcon} from "./SmallUIComponents";

const PlaygroundStepProductPicker = ({ products, setSmiles, reactionName }) => {
    useEffect(() => {
        // If there is only one product, there is no choice the user has to make
        if (products.length === 1) {
            setSmiles(products[0].smiles);
        }
    }, [products, setSmiles]);

    return (
        <>
            <ArrowWithReactionName reactionName={reactionName} />

            {/* list of clickable products */}
            <div className="mol-row">
                {products.length > 1 && products.map((product, index) => (
                    <React.Fragment key={product.smiles}>
                        <MolImage smiles={product.smiles} encoding={product.encoding} setSmiles={setSmiles} />
                        {index < products.length - 1 && {plusIcon}}
                    </React.Fragment>
                ))}
            </div>
            <p className="multi-product-text">
                This produces {products.length} products. Click on the one you want to analyze next.
            </p>
        </>
    );
}

export default PlaygroundStepProductPicker;