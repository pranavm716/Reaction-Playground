const PlaygroundStepReactionPicker = ({reactions}) => {
    return (
        <div className="reaction-picker">
            <p>Choose a reaction to run:</p>
            {reactions.map((reaction, index) => (
                <button key={index}>{reaction.name}</button>
            ))}
        </div>
    );
}
 
export default PlaygroundStepReactionPicker;