import { useState } from "react";
import ChemDraw from "../components/ChemDraw";

const Playground = () => {
    const [smiles, setSmiles] = useState("");
    const [loopStarted, setLoopStarted] = useState(false);

    const handleRunReactions = () => {
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
                    <ChemDraw setSmiles={setSmiles} />
                </div>
                <div className="action-panel">
                    {smiles ? 
                        <button className="run-reactions-button" onClick={handleRunReactions}>
                            Run Reactions
                        </button> 
                        : 
                        <p>Draw a molecule to get started!</p>
                    }
                </div>
            </div>
        </>
    );
}

export default Playground;