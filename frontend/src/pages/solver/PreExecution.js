import { useState } from "react";
import ChemDraw from "../../components/ChemDraw";
import MolImage from "../../components/MolImage";
import {
  CautionText,
  ClearSelectionButton,
  doubleArrowIcon,
} from "../../components/SmallUIComponents";

const noneSelectedText = <CautionText text="Not set" />;

const PreExecution = ({
  startingSmiles,
  setStartingSmiles,
  startingEncoding,
  targetSmiles,
  setTargetSmiles,
  targetEncoding,
  handleRunSolver,
}) => {
  const [preSolverSmiles, setPreSolverSmiles] = useState("");

  return (
    <>
      <p>
        This is the Solver. Here, you can draw a starting molecule and a target
        molecule and the solver will try to find a reaction path between them.
      </p>
      <div className="two-panel-content">
        <div>
          <ChemDraw setSmiles={setPreSolverSmiles} />
        </div>
        <div className="grid-container">
          {/* first row */}
          <p className="grid-starting-text">
            <b>Starting molecule</b>
          </p>

          <div className="grid-double-arrow-placeholder"></div>

          <p className="grid-target-text">
            <b>Target molecule</b>
          </p>

          {/* second row */}
          <div className="grid-starting-molecule">
            {startingEncoding ? (
              <MolImage smiles={startingSmiles} encoding={startingEncoding} />
            ) : (
              noneSelectedText
            )}
          </div>

          <div className="grid-double-arrow">{doubleArrowIcon}</div>

          <div className="grid-target-molecule">
            {targetEncoding ? (
              <MolImage smiles={targetSmiles} encoding={targetEncoding} />
            ) : (
              noneSelectedText
            )}
          </div>

          {/* third row */}
          <div className="grid-starting-button">
            {startingEncoding ? (
              <ClearSelectionButton onClick={() => setStartingSmiles("")} />
            ) : (
              preSolverSmiles &&
              !startingSmiles && (
                <button
                  onClick={() => setStartingSmiles(preSolverSmiles)}
                  className="primary-colored-button"
                >
                  Set as starting molecule
                </button>
              )
            )}
          </div>

          {startingSmiles && targetSmiles ? (
            <button
              onClick={handleRunSolver}
              className="primary-colored-button"
              style={{ height: "34px" }}
            >
              Find pathway
            </button>
          ) : (
            <div className="grid-double-arrow-placeholder"></div>
          )}

          <div className="grid-target-button">
            {targetEncoding ? (
              <ClearSelectionButton onClick={() => setTargetSmiles("")} />
            ) : (
              preSolverSmiles &&
              !targetSmiles && (
                <button
                  onClick={() => setTargetSmiles(preSolverSmiles)}
                  className="primary-colored-button"
                >
                  Set as target molecule
                </button>
              )
            )}
          </div>
        </div>
      </div>
      <p>Draw a starting and a target molecule to get started!</p>
    </>
  );
};

export default PreExecution;
