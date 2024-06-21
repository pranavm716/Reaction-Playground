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
    const [molImage, setMolImage] = useState(""); // base-64 encoded str of the current molecule
    const [reactions, setReactions] = useState([]); // list of valid reactions for the current molecule

    const handleStartLoop = async () => {
        await axios.get('/playground/step', { params: { smiles } })
        .then(res => {
            const [encoding, validReactions] = res.data;
            setMolImage(encoding);
            setReactions(validReactions);
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
                            : <PlaygroundStepReactionPicker reactions={reactions}/>
                    }
                </div>
            </div>
        </>
    );
}

export default Playground;