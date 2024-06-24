const MolImage = ({ smiles, encoding, setSmiles, isHighlighted }) => {
    // Supports both clickable and non clickable mol images
    const className = setSmiles ? "clickable-mol-image" : "mol-image";
    return (
        <img
            onClick={setSmiles ? () => setSmiles(smiles) : undefined}
            className={className}
            width={200}
            height={200}
            src={`data:image/png;base64,${encoding}`}
            alt={smiles}
            style={isHighlighted ? { boxShadow: '0 0 8px rgba(0, 0, 0, 0.3)' } : undefined}
        />
    );

    // TODO: Add context menu on right click
}

export default MolImage;