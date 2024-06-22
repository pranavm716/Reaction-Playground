const ArrowWithReactionName = ({reactionName}) => {
    return (
        <div style={{
            position: 'relative',
            display: 'flex',
            justifyContent: 'center',
            alignItems: 'center',
        }}>
            <img height={130} src={process.env.PUBLIC_URL + "down-arrow.png"}></img>
            <span style={{ width: 200, position: 'absolute', left: 60 }}>{reactionName}</span>
        </div>
    );
}
 
export default ArrowWithReactionName;