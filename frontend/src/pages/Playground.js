import { useState } from "react";
import ChemDraw from "../components/ChemDraw";
import MolImage from "../components/MolImage";

const Playground = () => {
    // before playground loop
    const [smiles, setSmiles] = useState("");
    const [loopStarted, setLoopStarted] = useState(false);

    // during playground loop
    const [molImage, setMolImage] = useState(""); // base-64 encoded str of the current molecule

    const handleStartLoop = () => {
        setLoopStarted(true);
    }

    return (
        <>
            <p>
                Welcome to Playground Mode! Here, you can draw a molecule and
                experiment with running common organic reactions on it.
            </p>
            <div className="two-panel-content">
                <div className="display-panel">
                    {loopStarted ? 
                        <MolImage smiles={smiles} encoding={molImage} />
                        : 
                        <ChemDraw setSmiles={setSmiles} />
                    }
                </div>
                <div className="action-panel">
                    {smiles && !loopStarted ? 
                        <button className="run-reactions-button" onClick={handleStartLoop}>
                            Run Reactions
                        </button> 
                        : !loopStarted ?
                            <p>Draw a molecule to get started!</p>
                            : null
                    }
                </div>
            </div>
        </>
    );
}

export default Playground;