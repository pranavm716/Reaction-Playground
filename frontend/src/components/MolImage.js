const MolImage = ({smiles, encoding}) => {
    // TODO: Add context menu
    return (
        <img className="mol-image" width={250} height={250} src={`data:image/png;base64,${encoding}`} alt={smiles} />
    );
}
 
export default MolImage;