const ReactionPicker = ({ reactions, setReactionPicked, molImage }) => {
    return (
        <>
            {molImage}
            {reactions.length ?
                <>
                    <p>Choose a reaction to run:</p>
                    {reactions.map((reaction, index) => (
                        <button className="reaction-button" key={index} onClick={() => setReactionPicked(reaction)}>
                            {reaction.name}
                        </button>
                    ))}
                </>
                :
                <p>There are no compatible reactions.</p>
            }
        </>
    );
}

export default ReactionPicker;