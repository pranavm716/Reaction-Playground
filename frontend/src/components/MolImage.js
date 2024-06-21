const MolImage = ({smiles, encoding}) => {
    // TODO: Add context menu
    return (
        <img src={`data:image/png;base64,${encoding}`} alt={smiles} />
    );
}
 
export default MolImage;