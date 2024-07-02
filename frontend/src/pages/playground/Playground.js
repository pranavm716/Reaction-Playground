import axios from "axios";
import { useEffect, useState } from "react";
import { useSearchParams } from "react-router-dom";
import ChemDraw from "../../components/ChemDraw";
import MolImage from "../../components/MolImage";
import { BackButton } from "../../components/SmallUIComponents";
import {
  MISSING_REACTANTS_PROMPTS_ENDPOINT,
  REACTION_MULTIPLE_REACTANTS_ENDPOINT,
  REACTION_SINGLE_REACTANT_ENDPOINT,
  START_ENDPOINT,
} from "../../endpoints";
import ExtraReactantPicker from "./ExtraReactantPicker";
import HistoryViewer from "./HistoryViewer";
import ProductPicker from "./ProductPicker";
import ReactionPicker from "./ReactionPicker";

const defaultMissingReactantMetadata = {
  prompts: null,
  smiles: null,
  encodings: null,
};

const defaultHistoryState = {
  curMolMetadata: { smiles: "", encoding: null },
  missingReactantMetadata: defaultMissingReactantMetadata,
  reactionPicked: null,
  productMetadata: null,
  productPicked: null,
};

const setCurHistoryAttribute = (setHistory, attribute, value) => {
  setHistory((prev) => [
    ...prev.slice(0, prev.length - 1),
    {
      ...prev[prev.length - 1],
      [attribute]: value,
    },
  ]);
};

// TODO: Figure out back button navigation for history mode
const Playground = () => {
  // before playground loop/ before a reaction is picked
  const [preLoopSmiles, setPreLoopSmiles] = useState(""); // smiles before explicitly starting the playground loop

  const [searchParams, setSearchParams] = useSearchParams({ smiles: "" });
  const smiles = searchParams.get("smiles"); // keeps track of the smiles for this loop
  const [history, setHistory] = useState([defaultHistoryState]); // history of reactions
  const curHistoryState = history[history.length - 1];

  const setNewSmiles = (productSmiles, productIndex) => {
    setProductPicked(productIndex);
    setSearchParams({ smiles: productSmiles }, { replace: true });
  };

  // metadata unique for step
  const [stepMetadata, setStepMetadata] = useState(null); // structure: {encoding -> mol img encoding, validReactions -> [list of valid reactions]}

  // functions to update the current molecule metadata in history
  const setCurMolMetadata = ({
    smiles = curHistoryState.curMolMetadata.smiles,
    encoding = curHistoryState.curMolMetadata.encoding,
  }) => {
    setCurHistoryAttribute(setHistory, "curMolMetadata", {
      smiles,
      encoding,
    });
  };

  const setReactionPicked = (reactionPicked) => {
    setCurHistoryAttribute(setHistory, "reactionPicked", reactionPicked);
  };

  const setMissingReactantMetadata = ({
    prompts = curHistoryState.missingReactantMetadata.prompts,
    smiles = curHistoryState.missingReactantMetadata.smiles,
    encodings = curHistoryState.missingReactantMetadata.encodings,
  }) => {
    setCurHistoryAttribute(setHistory, "missingReactantMetadata", {
      prompts,
      smiles,
      encodings,
    });
  };

  const setProductMetadata = (productMetadata) => {
    setCurHistoryAttribute(setHistory, "productMetadata", productMetadata);
  };

  const setProductPicked = (productPicked) => {
    setCurHistoryAttribute(setHistory, "productPicked", productPicked);
  };

  // start loop immediately if there are smiles from URL on initial page load
  useEffect(() => {
    const handleStepStart = async () => {
      await axios
        .get(START_ENDPOINT, { params: { smiles } })
        .then((res) => {
          const [encoding, validReactions] = res.data;
          setStepMetadata({ encoding, validReactions });

          setCurMolMetadata({ smiles, encoding });
        })
        .catch((error) => {
          alert(error.response.data.detail);
          setSearchParams({ smiles: "" });
        });
    };

    if (!smiles) {
      setStepMetadata(null);
      setHistory([defaultHistoryState]);
      return;
    }

    // append new entry to history, but only after the first loop
    if (curHistoryState.productMetadata) {
      setHistory((prev) => [...prev, defaultHistoryState]);
    }

    handleStepStart();

    // * This useEffect should only run when the current smiles changes
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [smiles]);

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
            setProductMetadata(res.data);
            if (products.length === 1) {
              setNewSmiles(products[0].smiles, 0);
            }
          });
      }
    };

    const reactionPicked = curHistoryState.reactionPicked;
    if (reactionPicked) handleStepReaction();

    // * This useEffect should only run when the current reactionPicked changes
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [curHistoryState.reactionPicked]);

  useEffect(() => {
    const handleStepReactionMultipleReactants = async () => {
      await axios
        .post(REACTION_MULTIPLE_REACTANTS_ENDPOINT, {
          smiles,
          extra_reactant_smiles: curHistoryState.missingReactantMetadata.smiles,
          reaction_key: curHistoryState.reactionPicked.reaction_key,
        })
        .then((res) => {
          const [extraReactantEncodings, products] = res.data;
          setMissingReactantMetadata({ encodings: extraReactantEncodings });
          setProductMetadata(products);

          if (products.length === 1) {
            setNewSmiles(products[0].smiles, 0);
          }
        })
        .catch((error) => {
          // provided reactants were invalid for this reaction
          setMissingReactantMetadata({
            smiles: defaultMissingReactantMetadata.smiles,
          });
          alert(error.response.data.detail);
        });
    };

    if (curHistoryState.missingReactantMetadata.smiles)
      handleStepReactionMultipleReactants();

    // * This useEffect should only run when the current missingReactantSmilesPicked changes
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [curHistoryState.missingReactantMetadata.smiles]);

  // If the user picks a reaction that requires multiple reactants but then decides to go back to the reaction picker
  const cancelMultipleReactants = () => {
    setReactionPicked(null);
    setMissingReactantMetadata(defaultMissingReactantMetadata);
  };

  // If the user wants to go back to the previous molecule, rewinding history
  const toPreviousMolecule = () => {
    setHistory((prev) => {
      setSearchParams(
        { smiles: prev[prev.length - 2].curMolMetadata.smiles },
        { replace: true },
      );
      return [...prev.slice(0, prev.length - 2), defaultHistoryState];
    });
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
              {curHistoryState.missingReactantMetadata.prompts &&
                !curHistoryState.missingReactantMetadata.encodings && (
                  <ExtraReactantPicker
                    molImage={molImage}
                    missingReactantPrompts={
                      curHistoryState.missingReactantMetadata.prompts
                    }
                    setMissingReactantSmilesPicked={(smiles) =>
                      setMissingReactantMetadata({ smiles })
                    }
                    reaction={curHistoryState.reactionPicked}
                    cancelMultipleReactants={cancelMultipleReactants}
                  />
                )}
              <>
                {curHistoryState.productMetadata ? (
                  <ProductPicker
                    products={curHistoryState.productMetadata}
                    setNewSmiles={setNewSmiles}
                    reaction={curHistoryState.reactionPicked}
                    molImage={molImage}
                    missingReactantSmilesPicked={
                      curHistoryState.missingReactantMetadata.smiles
                    }
                    missingReactantEncodings={
                      curHistoryState.missingReactantMetadata.encodings
                    }
                  />
                ) : (
                  !curHistoryState.missingReactantMetadata.prompts && (
                    <ReactionPicker
                      reactions={stepMetadata.validReactions}
                      setReactionPicked={setReactionPicked}
                      molImage={molImage}
                    />
                  )
                )}
                {history.length > 1 &&
                  !curHistoryState.missingReactantMetadata.prompts && (
                    <BackButton
                      text="Previous molecule"
                      onClick={toPreviousMolecule}
                    />
                  )}
              </>
            </div>
          ) : (
            <ChemDraw setSmiles={setPreLoopSmiles} />
          )}
        </div>
        <div className="history-panel">
          {stepMetadata ? (
            <>
              <p>
                <b>History</b>
              </p>
              <HistoryViewer history={history} />
            </>
          ) : preLoopSmiles ? (
            <div className="run-reactions-div">
              <button
                className="primary-colored-button"
                onClick={() => setSearchParams({ smiles: preLoopSmiles })}
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
