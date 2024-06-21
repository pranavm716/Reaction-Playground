import { useEffect } from "react";

const PlaygroundStepProductPicker = ({ products, setSmiles }) => {
    useEffect(() => {
        if (products.length === 1) {
            setSmiles(products[0].smiles);
        }
    }, [products, setSmiles]);

    return (
        // temp
        <div>
            {products.length > 1 && products.map(product => (
                <div key={product.smiles}>
                    <button onClick={() => setSmiles(product.smiles)}>{product.smiles}</button>
                </div>
            ))}
            <p>This produces {products.length} products. Click on the one you want to analyze next.</p>
        </div>
    );
}

export default PlaygroundStepProductPicker;