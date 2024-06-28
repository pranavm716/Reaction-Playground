import { useMolImageContextMenu } from "../hooks";

const MolImage = ({ smiles, encoding, onClick, isHighlighted }) => {
  // Supports both clickable and non clickable mol images
  const className = onClick ? "clickable-mol-image" : "mol-image";
  const { contextMenu, setIsOpen, setAnchorPoint } = useMolImageContextMenu(smiles);

  return (
    <>
      <img
        onClick={onClick || undefined}
        className={className}
        width={200}
        height={200}
        src={`data:image/png;base64,${encoding}`}
        alt={smiles}
        style={
          isHighlighted
            ? { boxShadow: "0 0 8px rgba(0, 0, 0, 0.3)" }
            : undefined
        }
        onContextMenu={(e) => {
          if (typeof document.hasFocus === "function" && !document.hasFocus())
            return;

          e.preventDefault();
          setAnchorPoint({ x: e.clientX, y: e.clientY });
          setIsOpen(true);
        }}
      />
      {contextMenu}
    </>
  );
};

export default MolImage;
