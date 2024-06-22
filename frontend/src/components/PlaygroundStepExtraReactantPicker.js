import ArrowWithReactionName from "./ArrowWithReactionName";
import ChemDraw from "./ChemDraw";
import PlusIcon from "./PlusIcon";
import React, { useState } from "react";
import Modal from 'react-modal';

const cautionIcon = <svg style={{
    paddingRight: '5px',
    position: 'relative',
    top: '5px',
}} xmlns="http://www.w3.org/2000/svg" height="24px" viewBox="0 -960 960 960" width="24px" fill="#f00">
    <path d="M109.23-160 480-800l370.77 640H109.23ZM178-200h604L480-720 178-200Zm302-55.38q10.46 0 17.54-7.08 7.08-7.08 7.08-17.54 0-10.46-7.08-17.54-7.08-7.08-17.54-7.08-10.46 0-17.54 7.08-7.08 7.08-7.08 17.54 0 10.46 7.08 17.54 7.08 7.08 17.54 7.08Zm-20-89.24h40v-200h-40v200ZM480-460Z" />
</svg>

const closeIcon = <svg style={{
    cursor: 'pointer',
}} xmlns="http://www.w3.org/2000/svg" height="20px" viewBox="0 -960 960 960" width="20px" fill="#000000">
    <path d="m252.62-217.23-35.39-35.39L444.62-480 217.23-707.38l35.39-35.39L480-515.38l227.38-227.39 35.39 35.39L515.38-480l227.39 227.38-35.39 35.39L480-444.62 252.62-217.23Z" />
</svg>

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
                        <PlusIcon />

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