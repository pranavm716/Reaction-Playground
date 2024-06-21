import { useEffect, useState } from "react";
import ChemDraw from "../components/ChemDraw";
import MolImage from "../components/MolImage";
import axios from 'axios';
import PlaygroundStepReactionPicker from "../components/PlaygroundStepReactionPicker";

const STEP_START_ENDPOINT = '/playground/step-start';
const STEP_REACTION_ENDPOINT = '/playground/step-reaction';

const Playground = () => {
    // before playground loop
    const [smiles, setSmiles] = useState("");

    // during playground loop
    const [stepMetadata, setStepMetadata] = useState(null); // metadata for the current step: {mol img encoding, list of valid reactions}
    const [reactionKeyPicked, setReactionKeyPicked] = useState(null); // key of the reaction picked by the user
    const [productsMetadata, setProductsMetadata] = useState(null); // metadata for the products of the current step: list of [{mol img encoding, smiles]}

    const handleStepStart = async () => {
        await axios.get(STEP_START_ENDPOINT, { params: { smiles } })
            .then(res => {
                const [encoding, validReactions] = res.data;
                setStepMetadata({ encoding, validReactions });
            })
            .catch(error => {
                console.log(error.response);
            })
    }

    const handleStepReaction = async () => {
        if (!reactionKeyPicked) return;
        await axios.get(STEP_REACTION_ENDPOINT, { params: { smiles, reaction_key: reactionKeyPicked } })
            .then(res => {
                setProductsMetadata(res.data);
            })
            .catch(error => {
                console.log(error.response);
            })
    }
    useEffect(() => { handleStepReaction() }, [reactionKeyPicked]);

    return (
        <>
            <p>
                Welcome to Playground Mode! Here, you can draw a molecule and
                experiment with running common organic reactions on it.
            </p>
            <div className="two-panel-content">
                <div className="action-panel">
                    {
                        productsMetadata ?
                            null
                            :
                            stepMetadata ?
                                <PlaygroundStepReactionPicker
                                    MolImageProp={<MolImage smiles={smiles} encoding={stepMetadata.encoding} />}
                                    reactions={stepMetadata.validReactions}
                                    setReactionKeyPicked={setReactionKeyPicked}
                                />
                                :
                                <ChemDraw setSmiles={setSmiles} />
                    }
                </div>
                <div className="history-panel">
                    {
                        stepMetadata ?
                            <p><b>History</b></p>
                            : smiles ?
                                <div className="run-reactions-div">
                                    <button onClick={handleStepStart}>
                                        Run Reactions
                                    </button>
                                </div>
                                :
                                <div className="run-reactions-div">
                                    <p>Draw a molecule to get started!</p>
                                </div>
                    }
                </div>
            </div>
        </>
    );
}

export default Playground;