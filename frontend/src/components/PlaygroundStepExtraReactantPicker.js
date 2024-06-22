import ArrowWithReactionName from "./ArrowWithReactionName";
import PlusIcon from "./PlusIcon";
import React from "react";

const PlaygroundStepExtraReactantPicker = ({ molImage, missingReactantPrompts, setMissingReactantSmilesPicked, reactionName }) => {
    return (
        <>
            <div className="mol-row">
                {molImage}
                {missingReactantPrompts.map(prompt => (
                    <React.Fragment key={prompt}>
                        <PlusIcon />

                        {/* TODO: Implement the chemdraw modal and output smiles */}
                        {/* TODO: Styles for this button */}
                        <button onClick={() => {}}>
                            {prompt}
                        </button>
                    </React.Fragment>
                ))}
            </div>

            <ArrowWithReactionName reactionName={reactionName} />
            . . .
        </>
    );
}

export default PlaygroundStepExtraReactantPicker;