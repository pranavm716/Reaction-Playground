import React, { useState } from "react";
import ChemDraw from "../components/ChemDraw";
import { cautionIcon, doubleArrowIcon } from "../components/SmallUIComponents";

const noneSelectedText = <div>
    {cautionIcon}
    Not set
</div>;

const Solver = () => {
    // before solver execution
    const [preLoopSmiles, setPreLoopSmiles] = useState('');

    const [startingSmiles, setStartingSmiles] = useState(null);
    const [startingEncoding, setStartingEncoding] = useState(null);

    const [targetSmiles, setTargetSmiles] = useState(null);
    const [targetEncoding, setTargetEncoding] = useState(null);

    // after solver execution

    return (
        <>
            <p>
                This is the Solver. Here, you can draw a starting molecule and a target molecule and
                the solver will try to find a reaction path between them.
            </p>
            <div className="two-panel-content">
                <div>
                    <ChemDraw setSmiles={setPreLoopSmiles} />
                </div>
                <div className="grid-container">
                    {/* first row */}
                    <p className="grid-starting-text"><b>Starting molecule</b></p>
                    <div className="grid-double-arrow-placeholder"></div>
                    <p className="grid-target-text"><b>Target molecule</b></p>

                    {/* second row */}
                    <div className="grid-starting-molecule">
                        {startingSmiles || noneSelectedText}
                    </div>
                    <div className="grid-double-arrow">{doubleArrowIcon}</div>
                    <div className="grid-target-molecule">
                        {targetSmiles || noneSelectedText}
                    </div>

                    {/* third row */}
                    <div className="grid-starting-button">
                        {preLoopSmiles && !startingSmiles &&
                            <button onClick={() => setStartingSmiles(preLoopSmiles)} className="primary-colored-button">
                                Set as starting molecule
                            </button>
                        }
                    </div>
                    {
                        startingSmiles && targetSmiles &&
                        <button className="primary-colored-button">
                            Find pathway
                        </button>
                    }
                    <div className="grid-button-placeholder"></div>
                    <div className="grid-target-button">
                        {preLoopSmiles && !targetSmiles &&
                            <button onClick={() => setTargetSmiles(preLoopSmiles)} className="primary-colored-button">
                                Set as target molecule
                            </button>
                        }
                    </div>
                </div>
            </div>
            <p>Draw 2 molecules to get started!</p>
        </>
    );
}

export default Solver;