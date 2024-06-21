import { useState } from "react";
import ChemDraw from "../components/ChemDraw";
import MolImage from "../components/MolImage";
import axios from 'axios';
import PlaygroundStepReactionPicker from "../components/PlaygroundStepReactionPicker";

const Playground = () => {
    // before playground loop
    const [smiles, setSmiles] = useState("");
    const [loopStarted, setLoopStarted] = useState(false);

    // during playground loop
    const [stepMetadata, setStepMetadata] = useState(null); // metadata for the current step (mol img encoding + valid reactions)


    const handleStartLoop = async () => {
        await axios.get('/playground/step', { params: { smiles } })
            .then(res => {
                const [encoding, validReactions] = res.data;
                setStepMetadata({ encoding, validReactions });
            })
            .catch(error => {
                console.log(error.response);
            })
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
                        <MolImage smiles={smiles} encoding={stepMetadata.encoding} />
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
                            : <PlaygroundStepReactionPicker reactions={stepMetadata.validReactions} />
                    }
                </div>
            </div>
        </>
    );
}

export default Playground;