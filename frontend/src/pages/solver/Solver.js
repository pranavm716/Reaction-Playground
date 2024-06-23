import { useState } from "react";
import PreExecution from "./PreExecution";

const Solver = () => {
    const [startingSmiles, setStartingSmiles] = useState(null);
    const [targetSmiles, setTargetSmiles] = useState(null);

    const handleRunSolver = async () => {
        
    }

    return (
        <PreExecution 
            startingSmiles={startingSmiles}
            targetSmiles={targetSmiles}
            setStartingSmiles={setStartingSmiles}
            setTargetSmiles={setTargetSmiles}
            handleRunSolver={handleRunSolver}
        />
    );
}

export default Solver;