import MolImage from "../../components/MolImage";

const SolverResults = ({
    startingSmiles,
    startingEncoding,
    targetSmiles,
    targetEncoding,
    solverResults,
}) => {
    return (
        <>
            <MolImage smiles={startingSmiles} encoding={startingEncoding} />
            <MolImage smiles={targetSmiles} encoding={targetEncoding} />
            <p>{solverResults.num_steps}</p>
        </>
    );
}

export default SolverResults;