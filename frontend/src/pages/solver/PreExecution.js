import React, { useState } from "react";
import ChemDraw from "../../components/ChemDraw";
import { cautionIcon, doubleArrowIcon } from "../../components/SmallUIComponents";
import axios from 'axios';
import MolImage from "../../components/MolImage";

const GET_MOL_IMAGE_ENDPOINT = '/solver/get-mol-image';
const noneSelectedText = <div>
    {cautionIcon}
    Not set
</div>;

const PreExecution = ({
    startingSmiles,
    setStartingSmiles,
    targetSmiles,
    setTargetSmiles,
    startingEncoding,
    setStartingEncoding,
    targetEncoding,
    setTargetEncoding,
    handleRunSolver,
}) => {
    const [preLoopSmiles, setPreLoopSmiles] = useState('');

    const handleSetSolverMolecule = async (isStartingMolecule) => {
        let setSmilesFn;
        let setEncodingFn;
        if (isStartingMolecule) {
            setSmilesFn = setStartingSmiles;
            setEncodingFn = setStartingEncoding;
        } else {
            setSmilesFn = setTargetSmiles;
            setEncodingFn = setTargetEncoding;
        }

        await axios.get(GET_MOL_IMAGE_ENDPOINT, { params: { smiles: preLoopSmiles } })
        .then(res => {
                setSmilesFn(preLoopSmiles);
                setEncodingFn(res.data);
            })
            .catch(error => {
                alert(error.response.data.detail);
                setSmilesFn("");
            })
    }

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
                        {
                            startingEncoding ?
                                <MolImage smiles={startingSmiles} encoding={startingEncoding} />
                                :
                                noneSelectedText
                        }
                    </div>
                    <div className="grid-double-arrow">{doubleArrowIcon}</div>
                    <div className="grid-target-molecule">
                        {
                            targetEncoding ?
                                <MolImage smiles={targetSmiles} encoding={targetEncoding} />
                                :
                                noneSelectedText
                        }
                    </div>

                    {/* third row */}
                    {/* TODO: add clear buttons for starting and target molecules */}
                    <div className="grid-starting-button">
                        {preLoopSmiles && !startingSmiles &&
                            <button onClick={() => handleSetSolverMolecule(true)} className="primary-colored-button">
                                Set as starting molecule
                            </button>
                        }
                    </div>
                    {
                        startingSmiles && targetSmiles &&
                        <button onClick={handleRunSolver} className="primary-colored-button">
                            Find pathway
                        </button>
                    }
                    <div className="grid-button-placeholder"></div>
                    <div className="grid-target-button">
                        {preLoopSmiles && !targetSmiles &&
                            <button onClick={() => handleSetSolverMolecule(false)} className="primary-colored-button">
                                Set as target molecule
                            </button>
                        }
                    </div>
                </div>
            </div>
            <p>Draw a starting and a target molecule to get started!</p>
        </>
    );
}

export default PreExecution;