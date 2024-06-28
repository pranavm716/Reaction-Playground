import ChemDraw from "../../components/ChemDraw";
import React, { useState } from "react";
import Modal from "react-modal";
import {
  CautionText,
  closeIcon,
  plusIcon,
  ArrowWithReactionInfo,
  BackButton,
} from "../../components/SmallUIComponents";

// TODO: Extract modal into a custom hook into a separate file (include other custom hooks too)
const ExtraReactantPicker = ({
  molImage,
  missingReactantPrompts,
  setMissingReactantSmilesPicked,
  reaction,
  cancelMultipleReactants,
}) => {
  const [modalIsOpen, setModalIsOpen] = useState(false);
  const [smiles, setSmiles] = useState(""); // Smiles of the current missing reactant picked

  // Modal styles and state management
  const customStyles = {
    content: {
      top: "50%",
      left: "50%",
      right: "auto",
      bottom: "auto",
      width: "60%",
      transform: "translate(-50%, -50%)",
    },
    closeIcon: {
      display: "flex",
      justifyContent: "flex-end",
    },
    header: {
      display: "flex",
      justifyContent: "space-between",
      alignItems: "center",
    },
    footer: {
      display: "flex",
      justifyContent: "center",
    },
  };

  const openModal = () => {
    setModalIsOpen(true);
  };

  const closeModal = () => {
    setModalIsOpen(false);
  };

  // Array and state management for the smiles of the missing reactants picked
  let missingReactantSmiles = new Array(missingReactantPrompts.length).fill(
    null
  );

  const handleUpdateMissingSmiles = (index) => {
    closeModal();
    missingReactantSmiles[index] = smiles;

    // If all missing reactants have been provided,
    // then update the state and let the parent's useEffect handle it from here
    if (missingReactantSmiles.every((smile) => smile !== null)) {
      setMissingReactantSmilesPicked(missingReactantSmiles);
    }
  };

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
            <Modal
              isOpen={modalIsOpen}
              onRequestClose={closeModal}
              contentLabel="Chem Draw Modal"
              style={customStyles}
            >
              <div style={customStyles.header}>
                {prompt} for {reaction.name}
                <div style={customStyles.closeIcon}>
                  {React.cloneElement(closeIcon, { onClick: closeModal })}
                </div>
              </div>

              <ChemDraw setSmiles={setSmiles} />

              <div style={customStyles.footer}>
                <button
                  onClick={() => {
                    handleUpdateMissingSmiles(index);
                  }}
                  className="primary-colored-button"
                >
                  Save
                </button>
              </div>
            </Modal>
          </React.Fragment>
        ))}
      </div>

      <ArrowWithReactionInfo
        reactionName={reaction.name}
        reactionDescriptionTooltip={reaction.description}
        tooltipId={`extra-reactants-${reaction.name}`}
      />
      <CautionText text="This reaction requires additional reactants." />
      <BackButton onClick={cancelMultipleReactants} text="Back to reactions" />
    </>
  );
};

export default ExtraReactantPicker;
