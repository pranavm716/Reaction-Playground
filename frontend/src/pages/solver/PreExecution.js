import React, { useEffect, useState } from "react";
import ChemDraw from "../../components/ChemDraw";
import { CautionText, doubleArrowIcon } from "../../components/SmallUIComponents";
import axios from 'axios';
import MolImage from "../../components/MolImage";

const GET_MOL_IMAGE_ENDPOINT = '/solver/get-mol-image';
const noneSelectedText = <CautionText text="Not set" />

const PreExecution = ({
    startingSmiles,
    setStartingSmiles,
    targetSmiles,
    setTargetSmiles,
    handleRunSolver,
}) => {
    const [preSolverSmiles, setPreSolverSmiles] = useState('');
    const [startingEncoding, setStartingEncoding] = useState(null);
    const [targetEncoding, setTargetEncoding] = useState(null);

    useEffect(() => {
        if (startingSmiles) {
            handleSetSolverMolecule(true);
        }
        if (targetSmiles) {
            handleSetSolverMolecule(false);
        }
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, []);

    const handleSetSolverMolecule = async (isStartingMolecule) => {
        let smiles;
        let setSmilesFn;
        let setEncodingFn;
        let matchesOtherMolecule = false;
        if (isStartingMolecule) {
            smiles = startingSmiles;
            setSmilesFn = setStartingSmiles;
            setEncodingFn = setStartingEncoding;
            matchesOtherMolecule = !startingSmiles && preSolverSmiles === targetSmiles;
        } else {
            smiles = targetSmiles;
            setSmilesFn = setTargetSmiles;
            setEncodingFn = setTargetEncoding;
            matchesOtherMolecule = !targetSmiles && preSolverSmiles === startingSmiles;
        }

        if (matchesOtherMolecule) {
            alert("Starting and target molecules cannot be the same!");
            return;
        }

        await axios.get(GET_MOL_IMAGE_ENDPOINT, { params: { smiles: smiles || preSolverSmiles } })
            .then(res => {
                if (!smiles) {
                    setSmilesFn(preSolverSmiles);
                }
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
                    <ChemDraw setSmiles={setPreSolverSmiles} />
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
                        {preSolverSmiles && !startingSmiles &&
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
                        {preSolverSmiles && !targetSmiles &&
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