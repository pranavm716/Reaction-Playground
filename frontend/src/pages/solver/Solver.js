import { useState, useEffect, useCallback } from "react";
import PreExecution from "./PreExecution";
import axios from 'axios';
import SolverResults from "./SolverResults";
import { useSearchParams } from 'react-router-dom';

const GET_MOL_IMAGE_ENDPOINT = '/solver/get-mol-image';
const RUN_SOLVER_ENDPOINT = '/solver/run';

const Solver = () => {
    // before execution
    const [searchParams, setSearchParams] = useSearchParams();

    const startingSmiles = searchParams.get('startingSmiles') || '';
    const [startingEncoding, setStartingEncoding] = useState(null);

    const targetSmiles = searchParams.get('targetSmiles') || '';
    const [targetEncoding, setTargetEncoding] = useState(null);

    // after execution
    const [solverResults, setSolverResults] = useState(null);

    // handles the display of the starting and target molecules, 
    // both from the url and from the chemdraw in the PreExecution component
    const handleSetSolverMolecule = useCallback(async (isStartingMolecule, preSolverSmiles) => {
        let smiles;
        let setSmilesFn;
        let setEncodingFn;
        let matchesOtherMolecule = false;
        if (isStartingMolecule) {
            smiles = startingSmiles;
            setSmilesFn = (smiles) => setSearchParams({ startingSmiles: smiles, targetSmiles });
            setEncodingFn = setStartingEncoding;
            matchesOtherMolecule = (startingSmiles && startingSmiles === targetSmiles) || preSolverSmiles === targetSmiles;
        } else {
            smiles = targetSmiles;
            setSmilesFn = (smiles) => setSearchParams({ startingSmiles, targetSmiles: smiles });
            setEncodingFn = setTargetEncoding;
            matchesOtherMolecule = (targetSmiles && targetSmiles === startingSmiles) || preSolverSmiles === startingSmiles;
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
    }, [setSearchParams, startingSmiles, targetSmiles])

    // Effect to update state when searchParams change
    useEffect(() => {
        if (startingSmiles === targetSmiles && startingSmiles !== '') {
            alert("Starting and target molecules cannot be the same!");
            setSearchParams({});
            return;
        }

        if (startingSmiles) {
            handleSetSolverMolecule(true);
        } else if (!startingSmiles) {
            setStartingEncoding(null);
        }

        if (targetSmiles) {
            handleSetSolverMolecule(false);
        } else if (!targetSmiles) {
            setTargetEncoding(null);
        }

        setSolverResults(null);
    }, [searchParams, setSearchParams, handleSetSolverMolecule, startingSmiles, targetSmiles, startingEncoding, targetEncoding]);

    // runs the solver and updates the solverResults state
    const handleRunSolver = async () => {
        await axios.get(RUN_SOLVER_ENDPOINT,
            {
                params: {
                    start_smiles: startingSmiles,
                    target_smiles: targetSmiles,
                }
            }
        )
            .then(res => {
                setSolverResults(res.data);
            })
    }

    if (solverResults) {
        return <SolverResults
            startingSmiles={startingSmiles}
            targetSmiles={targetSmiles}
            solverResults={solverResults}
        />;
    } else {
        return <PreExecution
            startingSmiles={startingSmiles}
            startingEncoding={startingEncoding}
            targetSmiles={targetSmiles}
            targetEncoding={targetEncoding}
            handleSetSolverMolecule={handleSetSolverMolecule}
            handleRunSolver={handleRunSolver}
        />;
    }
}

export default Solver;