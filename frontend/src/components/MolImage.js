import { MenuItem, SubMenu, ControlledMenu } from "@szhsin/react-menu";
import "@szhsin/react-menu/dist/index.css";
import "@szhsin/react-menu/dist/transitions/slide.css";
import { useState } from "react";
import { openExternalIcon, copyIcon } from "./SmallUIComponents";
import { useNavigate } from "react-router-dom";

const useReroute = (url, searchParams) => {
  const navigate = useNavigate();

  return () => {
    navigate(`${url}?${new URLSearchParams(searchParams).toString()}`);
  };
};

const useContextMenu = (smiles) => {
  const [isOpen, setIsOpen] = useState(false);
  const [anchorPoint, setAnchorPoint] = useState({ x: 0, y: 0 });

  // reroutes
  const navigateToPlayground = useReroute("/", { smiles: smiles });
  const navigateToSolverAsStarting = useReroute("/solver", {
    startingSmiles: smiles,
  });
  const navigateToSolverAsTarget = useReroute("/solver", {
    targetSmiles: smiles,
  });
  const navigateToClassifier = useReroute("/classifier", { smiles: smiles });

  const contextMenu = (
    <ControlledMenu
      anchorPoint={anchorPoint}
      state={isOpen ? "open" : "closed"}
      onClose={() => setIsOpen(false)}
      transition
    >
      <MenuItem onClick={() => navigator.clipboard.writeText(smiles)}>
        Copy SMILES
        {copyIcon}
      </MenuItem>
      <MenuItem onClick={navigateToPlayground}>
        Open in Playground
        {openExternalIcon}
      </MenuItem>
      <SubMenu label="Open in Solver">
        <MenuItem onClick={navigateToSolverAsStarting}>
          As starting molecule
          {openExternalIcon}
        </MenuItem>
        <MenuItem onClick={navigateToSolverAsTarget}>
          As target molecule
          {openExternalIcon}
        </MenuItem>
      </SubMenu>
      <MenuItem onClick={navigateToClassifier}>
        Open in Classifier
        {openExternalIcon}
      </MenuItem>
    </ControlledMenu>
  );

  return { contextMenu, setIsOpen, setAnchorPoint };
};

const MolImage = ({ smiles, encoding, onClick, isHighlighted }) => {
  // Supports both clickable and non clickable mol images
  const className = onClick ? "clickable-mol-image" : "mol-image";
  const { contextMenu, setIsOpen, setAnchorPoint } = useContextMenu(smiles);

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
