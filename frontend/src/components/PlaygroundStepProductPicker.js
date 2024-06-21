import { useEffect } from "react";

const PlaygroundStepProductPicker = ({ products, setSmiles, reactionName }) => {
    useEffect(() => {
        if (products.length === 1) {
            setSmiles(products[0].smiles);
        }
    }, [products, setSmiles]);

    return (
        <>
            <div style={{
                position: 'relative',
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
            }}>
                <img height={130} src={process.env.PUBLIC_URL + "down-arrow.png"}></img>
                <span style={{ width: 200, position: 'absolute', left: 60 }}>{reactionName}</span>
            </div>

            {products.length > 1 && products.map(product => (
                // TODO: change to image and generalize MolImage to handle clicks
                <div key={product.smiles}>
                    <button onClick={() => setSmiles(product.smiles)}>{product.smiles}</button>
                </div>
            ))}
            <p>This produces {products.length} products. Click on the one you want to analyze next.</p>
        </>
    );
}

export default PlaygroundStepProductPicker;