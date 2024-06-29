import axios from "axios";
import { useCallback, useEffect, useState } from "react";
import { useSearchParams } from "react-router-dom";
import ChemDraw from "../../components/ChemDraw";
import MolImage from "../../components/MolImage";
import {
  MISSING_REACTANTS_PROMPTS_ENDPOINT,
  REACTION_MULTIPLE_REACTANTS_ENDPOINT,
  REACTION_SINGLE_REACTANT_ENDPOINT,
  START_ENDPOINT,
} from "../../endpoints";
import ExtraReactantPicker from "./ExtraReactantPicker";
import ProductPicker from "./ProductPicker";
import ReactionPicker from "./ReactionPicker";

const defaultMissingReactantMetadata = {
  prompts: null,
  smiles: null,
  encodings: null,
};

// TODO: Add history and back button when there are no reactions (go back in history)
const Playground = () => {
  // before playground loop/ before a reaction is picked
  const [preLoopSmiles, setPreLoopSmiles] = useState(""); // smiles before explicitly starting the playground loop
  const [searchParams, setSearchParams] = useSearchParams({ smiles: "" });
  const smiles = searchParams.get("smiles"); // keeps track of the smiles for this loop
  const setSmiles = useCallback(
    (smiles) => setSearchParams({ smiles }),
    [setSearchParams],
  );

  // after a reaction is picked
  const [stepMetadata, setStepMetadata] = useState(null); // metadata for the current step: {encoding -> mol img encoding, validReactions -> [list of valid reactions]}
  const [reactionPicked, setReactionPicked] = useState(null); // reaction object picked by the user

  // if the picked reaction requires multiple reactants
  const [missingReactantMetadata, setMissingReactantMetadata] = useState(
    defaultMissingReactantMetadata,
  ); // metadata for the missing reactants of the current step: {prompts: [list of string prompts], smiles: [list of smiles], encodings: [list of encodings]}

  // if there are multiple products after the products are generated
  const [productsMetadata, setProductsMetadata] = useState(null); // metadata for the products of the current step: list of [{encoding: mol img encoding, smiles: smiles]}

  // start loop immediately if there are smiles from URL on initial page load
  useEffect(() => {
    const handleStepStart = async () => {
      await axios
        .get(START_ENDPOINT, { params: { smiles } })
        .then((res) => {
          const [encoding, validReactions] = res.data;
          setStepMetadata({ encoding, validReactions });
        })
        .catch((error) => {
          alert(error.response.data.detail);
          setSmiles("");
        });
    };

    if (!smiles) {
      setStepMetadata(null);
      return;
    }

    // clean up previous state so UI is rendered properly
    setMissingReactantMetadata(defaultMissingReactantMetadata);
    setProductsMetadata(null);

    handleStepStart();
  }, [smiles, setSmiles]);

  useEffect(() => {
    const handleStepReaction = async () => {
      if (reactionPicked.multiple_reactants_prompts) {
        // reaction has multiple reactants, need to prompt user for missing reactants
        await axios
          .get(MISSING_REACTANTS_PROMPTS_ENDPOINT, {
            params: { smiles, reaction_key: reactionPicked.reaction_key },
          })
          .then((res) => {
            setMissingReactantMetadata({ prompts: res.data });
          });
      } else {
        // reaction has only one reactant
        await axios
          .get(REACTION_SINGLE_REACTANT_ENDPOINT, {
            params: { smiles, reaction_key: reactionPicked.reaction_key },
          })
          .then((res) => {
            const products = res.data;
            if (products.length === 1) {
              setSmiles(products[0].smiles);
            } else {
              setProductsMetadata(res.data);
            }
          });
      }
    };

    if (reactionPicked) handleStepReaction();

    // * This useEffect should only run when reactionPicked changes, and not smiles
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [reactionPicked]);

  useEffect(() => {
    const handleStepReactionMultipleReactants = async () => {
      await axios
        .post(REACTION_MULTIPLE_REACTANTS_ENDPOINT, {
          smiles,
          extra_reactant_smiles: missingReactantMetadata.smiles,
          reaction_key: reactionPicked.reaction_key,
        })
        .then((res) => {
          // clean up previous state so UI is rendered properly
          setMissingReactantMetadata(defaultMissingReactantMetadata);

          const [extraReactantEncodings, products] = res.data;
          if (products.length === 1) {
            setSmiles(products[0].smiles);
          } else {
            setMissingReactantMetadata((prev) => ({
              ...prev,
              encodings: extraReactantEncodings,
            }));
            setProductsMetadata(products);
          }
        })
        .catch((error) => {
          // provided reactants were invalid for this reaction
          alert(error.response.data.detail);
        });
    };

    if (missingReactantMetadata.smiles) handleStepReactionMultipleReactants();

    // * This useEffect should only run when missingReactantSmilesPicked changes, and not smiles
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [missingReactantMetadata.smiles]);

  // If the user picks a reaction that requires multiple reactants but then decides to go back to the reaction picker
  const cancelMultipleReactants = () => {
    setReactionPicked(null);
    setMissingReactantMetadata(defaultMissingReactantMetadata);
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
              {missingReactantMetadata.prompts && (
                <ExtraReactantPicker
                  molImage={molImage}
                  missingReactantPrompts={missingReactantMetadata.prompts}
                  setMissingReactantSmilesPicked={(missingReactantSmiles) =>
                    setMissingReactantMetadata((prev) => ({
                      ...prev,
                      smiles: missingReactantSmiles,
                    }))
                  }
                  reaction={reactionPicked}
                  cancelMultipleReactants={cancelMultipleReactants}
                />
              )}
              {productsMetadata ? (
                <ProductPicker
                  products={productsMetadata}
                  setSmiles={setSmiles}
                  reaction={reactionPicked}
                  molImage={molImage}
                  missingReactantSmilesPicked={missingReactantMetadata.smiles}
                  missingReactantEncodings={missingReactantMetadata.encodings}
                />
              ) : (
                !missingReactantMetadata.prompts && (
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
                onClick={() => setSmiles(preLoopSmiles)}
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
