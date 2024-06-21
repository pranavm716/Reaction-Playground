const PlaygroundStepReactionPicker = ({ MolImageProp, reactions, setReactionKeyPicked }) => {
    return (
        <div className="reaction-picker">
            <p><b>Current molecule</b></p>
            {MolImageProp}

            {
                reactions.length ?
                    (
                        <>
                            <p>Choose a reaction to run:</p>
                            {reactions.map((reaction, index) => (
                                <button key={index} onClick={() => setReactionKeyPicked(reaction.reaction_key)}>
                                    {reaction.name}
                                </button>
                            ))}
                        </>
                    )
                    :
                    <p>There are no compatible reactions.</p>
            }

        </div>
    );
}

export default PlaygroundStepReactionPicker;