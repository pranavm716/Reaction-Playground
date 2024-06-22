import ArrowWithReactionName from "./ArrowWithReactionName";
import PlusIcon from "./PlusIcon";
import React from "react";

const cautionIcon = <svg style={{
    paddingRight: '5px',
    position: 'relative',
    top: '5px',
}} xmlns="http://www.w3.org/2000/svg" height="24px" viewBox="0 -960 960 960" width="24px" fill="#000000">
    <path d="M109.23-160 480-800l370.77 640H109.23ZM178-200h604L480-720 178-200Zm302-55.38q10.46 0 17.54-7.08 7.08-7.08 7.08-17.54 0-10.46-7.08-17.54-7.08-7.08-17.54-7.08-10.46 0-17.54 7.08-7.08 7.08-7.08 17.54 0 10.46 7.08 17.54 7.08 7.08 17.54 7.08Zm-20-89.24h40v-200h-40v200ZM480-460Z"/>
</svg>

const PlaygroundStepExtraReactantPicker = ({ molImage, missingReactantPrompts, setMissingReactantSmilesPicked, reactionName }) => {
    return (
        <>
            <div className="mol-row">
                {molImage}
                {missingReactantPrompts.map(prompt => (
                    <React.Fragment key={prompt}>
                        <PlusIcon />

                        {/* TODO: Implement the chemdraw modal and output smiles */}
                        <button className="extra-reactant-button" onClick={() => { }}>
                            {prompt}
                        </button>
                    </React.Fragment>
                ))}
            </div>

            <ArrowWithReactionName reactionName={reactionName} />
            <div>

                {cautionIcon}
                <span>
                This reaction requires additional reactants.

                </span>
            </div>
        </>
    );
}

export default PlaygroundStepExtraReactantPicker;