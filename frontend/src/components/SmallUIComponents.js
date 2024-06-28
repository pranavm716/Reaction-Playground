import React from "react";
import { Tooltip } from "react-tooltip";
import { Menu, MenuButton, MenuItem, SubMenu } from "@szhsin/react-menu";
import { useReroute } from "../hooks";
import "@szhsin/react-menu/dist/index.css";
import "@szhsin/react-menu/dist/transitions/slide.css";

export const closeIcon = (
  <svg
    style={{
      cursor: "pointer",
    }}
    xmlns="http://www.w3.org/2000/svg"
    height="20px"
    viewBox="0 -960 960 960"
    width="20px"
  >
    <path d="m252.62-217.23-35.39-35.39L444.62-480 217.23-707.38l35.39-35.39L480-515.38l227.38-227.39 35.39 35.39L515.38-480l227.39 227.38-35.39 35.39L480-444.62 252.62-217.23Z" />
  </svg>
);

export const plusIcon = (
  <svg
    style={{
      padding: "0 10px",
    }}
    xmlns="http://www.w3.org/2000/svg"
    height="24px"
    viewBox="0 -960 960 960"
    width="24px"
    fill="#000000"
  >
    <path d="M450-450H220v-60h230v-230h60v230h230v60H510v230h-60v-230Z" />
  </svg>
);

export const doubleArrowIcon = (
  <svg
    style={{
      paddingTop: "10px",
      margin: "0 15px",
    }}
    xmlns="http://www.w3.org/2000/svg"
    height="24px"
    viewBox="0 -960 960 960"
    width="24px"
    fill="#000000"
  >
    <path d="M383-480 200-664l56-56 240 240-240 240-56-56 183-184Zm264 0L464-664l56-56 240 240-240 240-56-56 183-184Z" />
  </svg>
);

export const openExternalIcon = (
  <svg
    style={{
      marginLeft: "15px",
      position: "relative",
    }}
    xmlns="http://www.w3.org/2000/svg"
    height="20px"
    viewBox="0 -960 960 960"
    width="24px"
    fill="#000000"
  >
    <path d="M212.31-140Q182-140 161-161q-21-21-21-51.31v-535.38Q140-778 161-799q21-21 51.31-21h252.3v60h-252.3q-4.62 0-8.46 3.85-3.85 3.84-3.85 8.46v535.38q0 4.62 3.85 8.46 3.84 3.85 8.46 3.85h535.38q4.62 0 8.46-3.85 3.85-3.84 3.85-8.46v-252.3h60v252.3Q820-182 799-161q-21 21-51.31 21H212.31Zm176.46-206.62-42.15-42.15L717.85-760H560v-60h260v260h-60v-157.85L388.77-346.62Z" />
  </svg>
);

export const copyIcon = (
  <svg
    style={{
      marginLeft: "15px",
      position: "relative",
    }}
    xmlns="http://www.w3.org/2000/svg"
    height="20px"
    viewBox="0 -960 960 960"
    width="20px"
    fill="#000000"
  >
    <path d="M362.31-260q-27.01 0-45.66-18.65Q298-297.3 298-324.31v-455.38q0-27.01 18.65-45.66Q335.3-844 362.31-844h359.38q27.01 0 45.66 18.65Q786-806.7 786-779.69v455.38q0 27.01-18.65 45.66Q748.7-260 721.69-260H362.31Zm0-52h359.38q4.62 0 8.46-3.85 3.85-3.84 3.85-8.46v-455.38q0-4.62-3.85-8.46-3.84-3.85-8.46-3.85H362.31q-4.62 0-8.46 3.85-3.85 3.84-3.85 8.46v455.38q0 4.62 3.85 8.46 3.84 3.85 8.46 3.85Zm-124 176q-27.01 0-45.66-18.65Q174-173.3 174-200.31v-507.38h52v507.38q0 4.62 3.85 8.46 3.84 3.85 8.46 3.85h411.38v52H238.31ZM350-312v-480 480Z" />
  </svg>
);

export const CautionText = ({ text }) => {
  return (
    <div>
      <svg
        style={{
          paddingRight: "5px",
          position: "relative",
          top: "5px",
        }}
        xmlns="http://www.w3.org/2000/svg"
        height="24px"
        viewBox="0 -960 960 960"
        width="24px"
        fill="#B00020"
      >
        <path d="M109.23-160 480-800l370.77 640H109.23ZM178-200h604L480-720 178-200Zm302-55.38q10.46 0 17.54-7.08 7.08-7.08 7.08-17.54 0-10.46-7.08-17.54-7.08-7.08-17.54-7.08-10.46 0-17.54 7.08-7.08 7.08-7.08 17.54 0 10.46 7.08 17.54 7.08 7.08 17.54 7.08Zm-20-89.24h40v-200h-40v200ZM480-460Z" />
      </svg>
      {text}
    </div>
  );
};

export const ArrowWithReactionInfo = ({
  reactionName,
  stepNumber,
  reactionDescriptionTooltip,
  tooltipId,
}) => {
  return (
    <div
      style={{
        position: "relative",
        display: "flex",
        justifyContent: "center",
        alignItems: "center",
        margin: "20px",
      }}
    >
      {stepNumber && (
        <span
          style={{
            position: "absolute",
            right: "50%",
            marginRight: "50px",
            whiteSpace: "nowrap",
          }}
        >
          <b>Step {stepNumber}</b>
        </span>
      )}
      <img
        height={130}
        src={process.env.PUBLIC_URL + "down-arrow.png"}
        alt="downward_arrow"
      ></img>
      <span
        data-tooltip-id={tooltipId}
        style={{
          position: "absolute",
          left: "50%",
          marginLeft: "50px",
          width: "200px",
        }}
      >
        {reactionName}
      </span>
      <Tooltip
        id={tooltipId}
        place="bottom"
        offset={12}
        content={reactionDescriptionTooltip}
        style={{ maxWidth: "350px", zIndex: 1000 }}
      />
    </div>
  );
};

export const BackButton = ({ text, onClick }) => {
  return (
    <button
      onClick={onClick}
      style={{
        height: "34px",
        marginTop: "20px",
        display: "inline-flex",
        alignItems: "center",
        justifyContent: "center",
        gap: "5px",
      }}
      className="primary-colored-button"
    >
      <svg
        xmlns="http://www.w3.org/2000/svg"
        height="23px"
        viewBox="0 -960 960 960"
        width="24px"
        fill="#fff"
      >
        <path d="m294.92-450 227.85 227.85L480-180 180-480l300-300 42.77 42.15L294.92-510H780v60H294.92Z" />
      </svg>
      {text}
    </button>
  );
};

export const ClearSelectionButton = ({ onClick }) => {
  return (
    <button onClick={onClick} className="clear-selection-button">
      {closeIcon}
      Clear selection
    </button>
  );
};

export const MolImageMenuItems = ({ smiles, onClassifier }) => {
  // reroutes
  const navigateToPlayground = useReroute("/", { smiles: smiles });
  const navigateToSolverAsStarting = useReroute("/solver", {
    startingSmiles: smiles,
  });
  const navigateToSolverAsTarget = useReroute("/solver", {
    targetSmiles: smiles,
  });
  const navigateToClassifier = useReroute("/classifier", { smiles: smiles });

  return (
    <>
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
      {!onClassifier && (
        <MenuItem onClick={navigateToClassifier}>
          Open in Classifier
          {openExternalIcon}
        </MenuItem>
      )}
    </>
  );
};

export const MolImageMenu = ({ smiles }) => {
  return (
    <Menu
      menuButton={
        <MenuButton className="primary-colored-button">Options</MenuButton>
      }
      transition
    >
      <MolImageMenuItems smiles={smiles} onClassifier={true} />
    </Menu>
  );
};
