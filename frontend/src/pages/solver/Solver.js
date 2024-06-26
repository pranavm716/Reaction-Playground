import { useState } from "react";
import PreExecution, { } from "./PreExecution";
import axios from 'axios';
import SolverResults from "./SolverResults";
import { useSearchParams } from 'react-router-dom';

const RUN_SOLVER_ENDPOINT = '/solver/run';

const Solver = () => {
    // before execution
    const [searchParams, setSearchParams] = useSearchParams({ startingSmiles: '', targetSmiles: '' });
    const startingSmiles = searchParams.get('startingSmiles');
    const targetSmiles = searchParams.get('targetSmiles');

    // after execution
    const [solverResults, setSolverResults] = useState(null);

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
            setStartingSmiles={(smiles) => setSearchParams({ startingSmiles: smiles, targetSmiles })}
            targetSmiles={targetSmiles}
            setTargetSmiles={(smiles) => setSearchParams({ startingSmiles, targetSmiles: smiles })}
            handleRunSolver={handleRunSolver}
        />;
    }
}

export default Solver;