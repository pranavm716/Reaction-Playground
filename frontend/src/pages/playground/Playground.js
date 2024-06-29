import { useCallback, useEffect, useState } from "react";
import ChemDraw from "../../components/ChemDraw";
import MolImage from "../../components/MolImage";
import axios from "axios";
import ReactionPicker from "./ReactionPicker";
import ProductPicker from "./ProductPicker";
import ExtraReactantPicker from "./ExtraReactantPicker";
import { useSearchParams } from "react-router-dom";

const START_ENDPOINT = "/playground/start";
const REACTION_SINGLE_REACTANT_ENDPOINT =
  "/playground/reaction/single-reactant";
const MISSING_REACTANTS_PROMPTS_ENDPOINT = "/playground/missing-reactants";
const REACTION_MULTIPLE_REACTANTS_ENDPOINT =
  "/playground/reaction/multiple-reactants";

// TODO: Add history and back button when there are no reactions (go back in history)
// TODO: Revisit all the no deps useEffects
const Playground = () => {
  // before playground loop
  const [preLoopSmiles, setPreLoopSmiles] = useState(""); // smiles before explicitly starting the playground loop
  const [searchParams, setSearchParams] = useSearchParams({ smiles: "" });

  // before a reaction is picked
  const smiles = searchParams.get("smiles"); // keeps track of the smiles for this loop

  // start loop immediately if there are smiles from URL on initial page load
  useEffect(() => {
    if (smiles) {
      handleStepStart(smiles);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const handleLoopStart = () => {
    handleStepStart(preLoopSmiles);
  };

  // after a reaction is picked
  const [stepMetadata, setStepMetadata] = useState(null); // metadata for the current step: {encoding -> mol img encoding, validReactions -> [list of valid reactions]}
  const [reactionPicked, setReactionPicked] = useState(null); // reaction object picked by the user

  // if the picked reaction requires multiple reactants
  const [missingReactantPrompts, setMissingReactantPrompts] = useState(null); // prompts for missing reactants: [list of string prompts]
  const [missingReactantSmilesPicked, setMissingReactantSmilesPicked] =
    useState(null); // smiles of the missing reactant picked by the user: [list of smiles]
  const [missingReactantEncodings, setMissingReactantEncodings] =
    useState(null); // encodings of the missing reactants picked by the user: [list of encodings]

  // after the products are generated
  const [productsMetadata, setProductsMetadata] = useState(null); // metadata for the products of the current step: list of [{encoding: mol img encoding, smiles: smiles]}

  // TODO: Reconsider making this a useEffect
  const handleStepStart = async (curSmiles) => {
    // clean up previous state so UI is rendered properly
    setProductsMetadata(null);
    setMissingReactantSmilesPicked(null);
    setMissingReactantEncodings(null);

    await axios
      .get(START_ENDPOINT, { params: { smiles: curSmiles } })
      .then((res) => {
        const [encoding, validReactions] = res.data;
        setStepMetadata({ encoding, validReactions });
        setSearchParams({ smiles: curSmiles });
      })
      .catch((error) => {
        alert(error.response.data.detail);
      });
  };

  const handleStepReaction = useCallback(async () => {
    if (!reactionPicked) return;

    if (reactionPicked.multiple_reactants_prompts) {
      // reaction has multiple reactants, need to prompt user for missing reactants
      await axios
        .get(MISSING_REACTANTS_PROMPTS_ENDPOINT, {
          params: { smiles, reaction_key: reactionPicked.reaction_key },
        })
        .then((res) => {
          setMissingReactantPrompts(res.data);
        });
    } else {
      // reaction has only one reactant
      await axios
        .get(REACTION_SINGLE_REACTANT_ENDPOINT, {
          params: { smiles, reaction_key: reactionPicked.reaction_key },
        })
        .then((res) => {
          setProductsMetadata(res.data);
        });
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [reactionPicked]);
  useEffect(() => {
    handleStepReaction();
  }, [handleStepReaction]);

  const handleStepReactionMultipleReactants = useCallback(async () => {
    if (!missingReactantSmilesPicked) return;

    await axios
      .post(REACTION_MULTIPLE_REACTANTS_ENDPOINT, {
        smiles,
        extra_reactant_smiles: missingReactantSmilesPicked,
        reaction_key: reactionPicked.reaction_key,
      })
      .then((res) => {
        // clean up previous state so UI is rendered properly
        setMissingReactantPrompts(null);

        const [extraReactantEncodings, products] = res.data;
        setMissingReactantEncodings(extraReactantEncodings);
        setProductsMetadata(products);
      })
      .catch((error) => {
        // provided reactants were invalid for this reaction
        alert(error.response.data.detail);
      });
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [missingReactantSmilesPicked]);
  useEffect(() => {
    handleStepReactionMultipleReactants();
  }, [handleStepReactionMultipleReactants]);

  const cancelMultipleReactants = () => {
    // If the user picks a reaction that requires multiple reactants but then decides to go back to the reaction picker
    setReactionPicked(null);
    setMissingReactantPrompts(null);
  };

  let molImage = null;
  if (stepMetadata) {
    molImage = <MolImage smiles={smiles} encoding={stepMetadata.encoding} />;
  }

  return (
    <>
      {!stepMetadata && (
        <p>
          Welcome to Reaction Playground! In playground mode, you can draw a
          molecule and experiment with running common organic reactions on it.
        </p>
      )}
      <div className="two-panel-content">
        <div className="action-panel">
          {stepMetadata ? (
            <div className="reaction-or-product-picker">
              <p>
                <b>Current molecule</b>
              </p>
              {missingReactantPrompts && (
                <ExtraReactantPicker
                  molImage={molImage}
                  missingReactantPrompts={missingReactantPrompts}
                  setMissingReactantSmilesPicked={
                    setMissingReactantSmilesPicked
                  }
                  reaction={reactionPicked}
                  cancelMultipleReactants={cancelMultipleReactants}
                />
              )}
              {productsMetadata ? (
                <ProductPicker
                  products={productsMetadata}
                  handleStepStart={handleStepStart}
                  reaction={reactionPicked}
                  molImage={molImage}
                  missingReactantSmilesPicked={missingReactantSmilesPicked}
                  missingReactantEncodings={missingReactantEncodings}
                />
              ) : (
                !missingReactantPrompts && (
                  <ReactionPicker
                    reactions={stepMetadata.validReactions}
                    setReactionPicked={setReactionPicked}
                    molImage={molImage}
                  />
                )
              )}
            </div>
          ) : (
            <ChemDraw setSmiles={setPreLoopSmiles} />
          )}
        </div>
        <div className="history-panel">
          {stepMetadata ? (
            <p>
              <b>History</b>
            </p>
          ) : preLoopSmiles ? (
            <div className="run-reactions-div">
              <button
                className="primary-colored-button"
                onClick={handleLoopStart}
              >
                Run Reactions
              </button>
            </div>
          ) : (
            <div className="run-reactions-div">
              <p>Draw a molecule to get started!</p>
            </div>
          )}
        </div>
      </div>
    </>
  );
};

export default Playground;
