import MolImage from "../../components/MolImage";
import { doubleArrowIcon } from "../../components/SmallUIComponents";

const SolverResults = ({
    startingSmiles,
    startingEncoding,
    targetSmiles,
    targetEncoding,
    solverResults,
}) => {
    return (
        <>
            <div>
                <p>
                    The goal is to find a reaction pathway that converts
                    the starting molecule into the target molecule.
                </p>
                <div className="solver-initial-display">
                    <MolImage smiles={startingSmiles} encoding={startingEncoding} />
                    {doubleArrowIcon}
                    <MolImage smiles={targetSmiles} encoding={targetEncoding} />
                </div>
                {
                    solverResults.path_found ?
                        <p>
                            {solverResults.num_steps} step synthetic pathway found!
                        </p>
                        :
                        <p>
                            {/* TODO: make the max num solver steps and multi step react mode as frontend configs  */}
                            No synthetic pathway found in 15 steps. It is also possible that no synthetic
                            pathway exists at all.
                        </p>
                }
            </div>
            {
                solverResults.path_found &&
                <div className="solver-results-display">
                    <p>This is your starting molecule:</p>
                </div>
            }
        </>
    );
}

export default SolverResults;