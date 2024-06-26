import React, { useCallback, useEffect } from "react";
import MolImage from "../../components/MolImage";
import { ArrowWithReactionInfo, plusIcon } from "../../components/SmallUIComponents";

const ProductPicker = ({ products, handleStepStart, reactionName, molImage, missingReactantSmilesPicked, missingReactantEncodings }) => {
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

    const pickProduct = useCallback((index) => {
        handleStepStart(products[index].smiles);
    }, [products, handleStepStart]);

    // If there is only one product, there is no choice the user has to make
    useEffect(() => {
        if (products.length === 1) {
            pickProduct(0);
        }
    }, [products, pickProduct]);

    return (
        <>
            {reactantMolRow}
            <ArrowWithReactionInfo reactionName={reactionName} />

            {/* list of clickable products - only shows when there is more than 1 */}
            {
                products.length > 1 &&
                <>
                    <div className="mol-row">
                        {products.map((product, index) => (
                            <React.Fragment key={product.smiles}>
                                <MolImage
                                    smiles={product.smiles}
                                    encoding={product.encoding}
                                    onClick={() => pickProduct(index)}
                                />
                                {index < products.length - 1 && plusIcon}
                            </React.Fragment>
                        ))}
                    </div>
                    <p className="multi-product-text">
                        This produces {products.length} products. Click on the one you want to analyze next.
                    </p>
                </>
            }
        </>
    );
}

export default ProductPicker;