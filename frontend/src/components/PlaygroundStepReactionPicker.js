const PlaygroundStepReactionPicker = ({ reactions, setReactionKeyPicked }) => {
    return (
        <div className="reaction-picker">
            <p>Choose a reaction to run:</p>
            {reactions.map((reaction, index) => (
                <button key={index} onClick={() => setReactionKeyPicked(reaction.reaction_key)}>
                    {reaction.name}
                </button>
            ))}
        </div>
    );
}

export default PlaygroundStepReactionPicker;