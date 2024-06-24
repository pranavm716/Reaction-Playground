import React from "react";
import MolImage from "../../components/MolImage";
import { ArrowWithReactionInfo, doubleArrowIcon } from "../../components/SmallUIComponents";
import { plusIcon } from "../../components/SmallUIComponents";

const SolverResults = ({
    startingSmiles,
    startingEncoding,
    targetSmiles,
    targetEncoding,
    solverResults,
}) => {
    const startMol = <MolImage smiles={startingSmiles} encoding={startingEncoding} />;
    const targetMol = <MolImage smiles={targetSmiles} encoding={targetEncoding} />;

    return (
        <>
            <div>
                <p>
                    The goal is to find a reaction pathway that converts
                    the starting molecule into the target molecule.
                </p>
                <div className="solver-initial-display">
                    {startMol}
                    {doubleArrowIcon}
                    {targetMol}
                </div>
                {
                    solverResults.path_found ?
                        <p>
                            <b style={{ color: 'green' }}>
                                {solverResults.num_steps} step synthetic pathway found!
                            </b>
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
                <div style={{ marginTop: '25px' }}>
                    <p>This is your starting molecule:</p>
                    <div className="solver-results-display">
                        {startMol}
                        {
                            // map over the synthetic pathway steps
                            [...Array(solverResults.num_steps).keys()].map(step => {
                                return (
                                    <div key={step}>
                                        <ArrowWithReactionInfo reactionName={solverResults.reaction_names[step]} stepNumber={step + 1} />
                                        <div style={{ display: 'flex', alignItems: 'center', position: 'relative' }}>
                                            <div className="solver-mol-row">

                                                {/* map over the products for this step */}
                                                {solverResults.solver_image_metadata[step].map((product, index) => (
                                                    <React.Fragment key={product.smiles}>
                                                        <MolImage
                                                            smiles={solverResults.solver_image_metadata[step][index].smiles}
                                                            encoding={solverResults.solver_image_metadata[step][index].encoding}
                                                            isHighlighted={solverResults.solver_image_metadata[step].length > 1 && solverResults.choice_pathway[step] === index}
                                                        />
                                                        {index < solverResults.solver_image_metadata[step].length - 1 && plusIcon}
                                                    </React.Fragment>
                                                ))}
                                            </div>

                                            {
                                                solverResults.solver_image_metadata[step].length > 1 &&
                                                <p style={{ marginLeft: '40px', position: 'absolute', left: '100%', width: '220px' }}>
                                                    This produces {solverResults.solver_image_metadata[step].length} products.
                                                    Pick product #{solverResults.choice_pathway[step] + 1}.
                                                </p>
                                            }
                                        </div>
                                    </div>
                                );
                            })
                        }
                    </div>
                    <p>This is your target molecule.</p>
                </div>
            }
            <button onClick={() => window.location.reload()} style={{ marginTop: '20px' }} className="primary-colored-button">Back to Solver</button>
        </>
    );
}

export default SolverResults;