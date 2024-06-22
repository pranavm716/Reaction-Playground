import { useEffect, useState } from "react";
import ChemDraw from "../components/ChemDraw";
import MolImage from "../components/MolImage";
import axios from 'axios';
import PlaygroundStepReactionPicker from "../components/PlaygroundStepReactionPicker";
import PlaygroundStepProductPicker from "../components/PlaygroundStepProductPicker";
import PlaygroundStepExtraReactantPicker from "../components/PlaygroundStepExtraReactantPicker";

const START_ENDPOINT = '/playground/start';
const REACTION_SINGLE_REACTANT_ENDPOINT = '/playground/reaction/single-reactant';
const MISSING_REACTANTS_PROMPTS_ENDPOINT = '/playground/missing-reactants';
const REACTION_MULTIPLE_REACTANTS_ENDPOINT = '/playground/reaction/multiple-reactants';

const Playground = () => {
    // before playground loop
    const [preLoopSmiles, setPreLoopSmiles] = useState(""); // smiles before explicitly starting the playground loop

    const handleLoopStart = () => {
        setSmiles(preLoopSmiles); // useEffect for smiles will handle the loop from here
    }

    // during playground loop
    const [smiles, setSmiles] = useState("");

    const [stepMetadata, setStepMetadata] = useState(null); // metadata for the current step: {encoding -> mol img encoding, validReactions -> [list of valid reactions]}
    const [reactionPicked, setReactionPicked] = useState(null); // reaction object picked by the user

    const [productsMetadata, setProductsMetadata] = useState(null); // metadata for the products of the current step: list of [{encoding: mol img encoding, smiles: smiles]}

    const [missingReactantPrompts, setMissingReactantPrompts] = useState(null); // prompts for missing reactants: [list of string prompts]
    const [missingReactantSmilesPicked, setMissingReactantSmilesPicked] = useState(null); // smiles of the missing reactant picked by the user: [list of smiles]
    const [missingReactantEncodings, setMissingReactantEncodings] = useState(null); // encodings of the missing reactants picked by the user: [list of encodings]

    const handleStepStart = async () => {
        if (!smiles) return;

        // clean up previous state so UI is rendered properly
        setProductsMetadata(null);
        setMissingReactantSmilesPicked(null);
        setMissingReactantEncodings(null);

        await axios.get(START_ENDPOINT, { params: { smiles } })
            .then(res => {
                const [encoding, validReactions] = res.data;
                setStepMetadata({ encoding, validReactions });
            })
            .catch(error => {
                alert(error.response.data.detail);
                setSmiles("");
            })
    }
    useEffect(() => { handleStepStart() }, [smiles]);

    const handleStepReaction = async () => {
        if (!reactionPicked) return;

        if (reactionPicked.multiple_reactants_prompts) {
            // reaction has multiple reactants
            await axios.get(MISSING_REACTANTS_PROMPTS_ENDPOINT, { params: { smiles, reaction_key: reactionPicked.reaction_key } })
                .then(res => {
                    setMissingReactantPrompts(res.data);
                })
        } else {
            // reaction has only one reactant
            await axios.get(REACTION_SINGLE_REACTANT_ENDPOINT, { params: { smiles, reaction_key: reactionPicked.reaction_key } })
                .then(res => {
                    setProductsMetadata(res.data.map(product => ({ encoding: product[0], smiles: product[1] })));
                })
        }
    }
    useEffect(() => { handleStepReaction() }, [reactionPicked]);

    const handleStepReactionMultipleReactants = async () => {
        if (!missingReactantSmilesPicked) return;

        await axios.post(REACTION_MULTIPLE_REACTANTS_ENDPOINT,
            {
                smiles: smiles,
                extra_reactant_smiles: missingReactantSmilesPicked,
                reaction_key: reactionPicked.reaction_key
            }
        )
            .then(res => {
                // clean up previous state so UI is rendered properly
                setMissingReactantPrompts(null);

                const [extraReactantEncodings, products] = res.data;
                setMissingReactantEncodings(extraReactantEncodings);
                setProductsMetadata(products.map(product => ({ encoding: product[0], smiles: product[1] })));
            })
            .catch(error => {
                // provided reactants were invalid for this reaction
                alert(error.response.data.detail);
            })
    }
    useEffect(() => { handleStepReactionMultipleReactants() }, [missingReactantSmilesPicked]);

    let molImage = null;
    if (stepMetadata) {
        molImage = <MolImage smiles={smiles} encoding={stepMetadata.encoding} />
    }

    return (
        <>
            {
                !stepMetadata && <p>
                    Welcome to Reaction Playground! In playground mode, you can draw a molecule and
                    experiment with running common organic reactions on it.
                </p>
            }
            <div className="two-panel-content">
                <div className="action-panel">
                    {
                        stepMetadata ?
                            <div className="reaction-or-product-picker">
                                <p><b>Current molecule</b></p>
                                {
                                    missingReactantPrompts &&
                                    <PlaygroundStepExtraReactantPicker
                                        molImage={molImage}
                                        missingReactantPrompts={missingReactantPrompts}
                                        setMissingReactantSmilesPicked={setMissingReactantSmilesPicked}
                                        reactionName={reactionPicked.name}
                                    />
                                }
                                {
                                    productsMetadata ?
                                        <PlaygroundStepProductPicker
                                            products={productsMetadata}
                                            setSmiles={setSmiles}
                                            reactionName={reactionPicked.name}
                                            molImage={molImage}
                                            missingReactantSmilesPicked={missingReactantSmilesPicked}
                                            missingReactantEncodings={missingReactantEncodings}
                                        />
                                        : !missingReactantPrompts ?
                                            <PlaygroundStepReactionPicker
                                                reactions={stepMetadata.validReactions}
                                                setReactionPicked={setReactionPicked}
                                                molImage={molImage}
                                            /> : null
                                }
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
                                    <button className="primary-colored-button" onClick={handleLoopStart}>
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