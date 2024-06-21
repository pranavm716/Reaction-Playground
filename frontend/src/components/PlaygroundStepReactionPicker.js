const PlaygroundStepReactionPicker = ({ reactions, setReactionKeyPicked }) => {
    return (
        <>
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
        </>
    );
}

export default PlaygroundStepReactionPicker;