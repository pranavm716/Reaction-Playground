import React from "react";
import {
  ArrowWithReactionInfo,
  BackButton,
  CautionText,
  plusIcon,
} from "../../components/SmallUIComponents";
import { useExtraReactantModal } from "../../hooks";

const ExtraReactantPicker = ({
  molImage,
  missingReactantPrompts,
  setMissingReactantSmilesPicked,
  reaction,
  cancelMultipleReactants,
}) => {
  // Array and state management for the smiles of the missing reactants picked
  let missingReactantSmiles = new Array(missingReactantPrompts.length).fill(
    null,
  );

  const { openModal, modals } = useExtraReactantModal(
    missingReactantSmiles,
    missingReactantPrompts,
    reaction.name,
    setMissingReactantSmilesPicked,
  );

  return (
    <>
      <div className="mol-row">
        {molImage}
        {missingReactantPrompts.map((prompt, index) => (
          <React.Fragment key={prompt}>
            {plusIcon}

            <button className="extra-reactant-button" onClick={openModal}>
              {prompt}
            </button>

            {modals[index]}
          </React.Fragment>
        ))}
      </div>

      <ArrowWithReactionInfo
        reactionName={reaction.name}
        reactionDescriptionTooltip={reaction.description}
      />
      <CautionText text="This reaction requires additional reactants." />
      <BackButton onClick={cancelMultipleReactants} text="Back to reactions" />
    </>
  );
};

export default ExtraReactantPicker;
