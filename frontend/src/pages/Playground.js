import { useEffect, useState } from "react";
import ChemDraw from "../components/ChemDraw";
import MolImage from "../components/MolImage";
import axios from 'axios';
import PlaygroundStepReactionPicker from "../components/PlaygroundStepReactionPicker";
import PlaygroundStepProductPicker from "../components/PlaygroundStepProductPicker";

const STEP_START_ENDPOINT = '/playground/step-start';
const STEP_REACTION_ENDPOINT = '/playground/step-reaction';

const Playground = () => {
    // before playground loop
    const [preLoopSmiles, setPreLoopSmiles] = useState("COC(C)=O"); // smiles before explicitly starting the playground loop

    const handleLoopStart = () => {
        setSmiles(preLoopSmiles); // useEffect for smiles will handle the loop from here
    }

    // during playground loop
    const [smiles, setSmiles] = useState("");
    const [stepMetadata, setStepMetadata] = useState(null); // metadata for the current step: {encoding -> mol img encoding, validReactions -> list of valid reactions}
    const [reactionPicked, setReactionPicked] = useState(null); // reaction object picked by the user
    const [productsMetadata, setProductsMetadata] = useState(null); // metadata for the products of the current step: list of [{encoding: mol img encoding, smiles -> smiles]}

    const handleStepStart = async () => {
        if (!smiles) return;

        // clean up previous state so UI is rendered properly
        setProductsMetadata(null);

        await axios.get(STEP_START_ENDPOINT, { params: { smiles } })
            .then(res => {
                const [encoding, validReactions] = res.data;
                setStepMetadata({ encoding, validReactions });
            })
            .catch(error => {
                console.log(error.response);
            })
    }
    useEffect(() => { handleStepStart() }, [smiles]);

    const handleStepReaction = async () => {
        if (!reactionPicked) return;
        await axios.get(STEP_REACTION_ENDPOINT, { params: { smiles, reaction_key: reactionPicked.reaction_key } })
            .then(res => {
                setProductsMetadata(res.data.map(product => ({ encoding: product[0], smiles: product[1] })));
            })
            .catch(error => {
                console.log(error.response);
            })
    }
    useEffect(() => { handleStepReaction() }, [reactionPicked]);

    let molImage = null;
    if (stepMetadata) {
        molImage = <>
            <p><b>Current molecule</b></p>
            <MolImage smiles={smiles} encoding={stepMetadata.encoding} />
        </>
    }

    return (
        <>
            {
                !stepMetadata && <p>
                    Welcome to Playground Mode! Here, you can draw a molecule and
                    experiment with running common organic reactions on it.
                </p>
            }
            <div className="two-panel-content">
                <div className="action-panel">
                    {
                        productsMetadata ?
                            <div className="reaction-or-product-picker">
                                {molImage}
                                <PlaygroundStepProductPicker 
                                    products={productsMetadata}
                                    setSmiles={setSmiles} 
                                    reactionName={reactionPicked.name}
                                />
                            </div>
                            :
                            stepMetadata ?
                                <div className="reaction-or-product-picker">
                                    {molImage}
                                    <PlaygroundStepReactionPicker
                                        reactions={stepMetadata.validReactions}
                                        setReactionPicked={setReactionPicked}
                                    />
                                </div>
                                :
                                <ChemDraw setSmiles={setPreLoopSmiles} />
                    }
                </div>
                <div className="history-panel">
                    {
                        stepMetadata ?
                            <p><b>History</b></p>
                            : preLoopSmiles ?
                                <div className="run-reactions-div">
                                    <button onClick={handleLoopStart}>
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