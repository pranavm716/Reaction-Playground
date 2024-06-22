import ChemDraw from "./ChemDraw";
import React, { useState } from "react";
import Modal from 'react-modal';
import { cautionIcon, closeIcon, plusIcon, ArrowWithReactionName } from "./SmallUIComponents";

const PlaygroundStepExtraReactantPicker = ({ molImage, missingReactantPrompts, setMissingReactantSmilesPicked, reactionName }) => {
    const [modalIsOpen, setModalIsOpen] = useState(false);
    const customStyles = {
        content: {
            top: '50%',
            left: '50%',
            right: 'auto',
            bottom: 'auto',
            width: '60%',
            transform: 'translate(-50%, -50%)',

        },
        closeIcon: {
            display: 'flex',
            justifyContent: 'flex-end',
        },
        header: {
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
        }
    };
    const openModal = () => {
        setModalIsOpen(true);
    }

    const closeModal = () => {
        setModalIsOpen(false);
    }

    return (
        <>
            <div className="mol-row">
                {molImage}
                {missingReactantPrompts.map(prompt => (
                    <React.Fragment key={prompt}>
                        {plusIcon}

                        {/* TODO: Output smiles */}
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
                                {prompt} for {reactionName}
                                <div style={customStyles.closeIcon}>
                                    {React.cloneElement(closeIcon, { onClick: closeModal })}
                                </div>
                            </div>
                            <ChemDraw setSmiles={setMissingReactantSmilesPicked} />
                        </Modal>
                    </React.Fragment>
                ))}
            </div>

            <ArrowWithReactionName reactionName={reactionName} />
            <div>
                {cautionIcon}
                This reaction requires additional reactants.
            </div>
        </>
    );
}

export default PlaygroundStepExtraReactantPicker;