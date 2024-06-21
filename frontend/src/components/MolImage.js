const MolImage = ({smiles, encoding, setSmiles}) => {
    // Supports both clickable and non clickable mol images
    const className = setSmiles ? "clickable-mol-image" : "mol-image";
    return (
        <img onClick={() => setSmiles(smiles)} className={className} width={200} height={200} src={`data:image/png;base64,${encoding}`} alt={smiles} />
    );
    
    // TODO: Add context menu on right click
}
 
export default MolImage;