const PlaygroundStepReactionPicker = ({ reactions, setReactionPicked }) => {
    return (
        <>
            {
                reactions.length ?
                (
                    <>
                        <p>Choose a reaction to run:</p>
                        {reactions.map((reaction, index) => (
                            <button key={index} onClick={() => setReactionPicked(reaction)}>
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