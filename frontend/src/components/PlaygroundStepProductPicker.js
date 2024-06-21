import React, { useEffect } from "react";
import MolImage from "./MolImage";

const plusIcon = <svg style={{
    padding: "0 10px",
}} xmlns="http://www.w3.org/2000/svg" height="24px" viewBox="0 -960 960 960" width="24px" fill="#000000">
    <path d="M450-450H220v-60h230v-230h60v230h230v60H510v230h-60v-230Z" />
</svg>


const PlaygroundStepProductPicker = ({ products, setSmiles, reactionName }) => {
    useEffect(() => {
        // If there is only one product, there is no choice the user has to make
        if (products.length === 1) {
            setSmiles(products[0].smiles);
        }
    }, [products, setSmiles]);

    return (
        <>
            {/* down arrow + reaction name */}
            <div style={{
                position: 'relative',
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
            }}>
                <img height={130} src={process.env.PUBLIC_URL + "down-arrow.png"}></img>
                <span style={{ width: 200, position: 'absolute', left: 60 }}>{reactionName}</span>
            </div>

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

export default PlaygroundStepProductPicker;