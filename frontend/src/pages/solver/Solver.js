import { useState } from "react";
import PreExecution from "./PreExecution";
import axios from 'axios';
import SolverResults from "./SolverResults";

const RUN_SOLVER_ENDPOINT = '/solver/run';

const Solver = () => {
    // before execution
    const [startingSmiles, setStartingSmiles] = useState(null);
    const [targetSmiles, setTargetSmiles] = useState(null);

    const [startingEncoding, setStartingEncoding] = useState(null);
    const [targetEncoding, setTargetEncoding] = useState(null);

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
            startingEncoding={startingEncoding}
            targetSmiles={targetSmiles}
            targetEncoding={targetEncoding}
            solverResults={solverResults}
        />
    } else {
        return <PreExecution
            startingSmiles={startingSmiles}
            setStartingSmiles={setStartingSmiles}
            targetSmiles={targetSmiles}
            setTargetSmiles={setTargetSmiles}
            startingEncoding={startingEncoding}
            setStartingEncoding={setStartingEncoding}
            targetEncoding={targetEncoding}
            setTargetEncoding={setTargetEncoding}
            handleRunSolver={handleRunSolver}
        />;
    }
}

export default Solver;