import { useState, useEffect, useCallback } from "react";
import PreExecution from "./PreExecution";
import axios from "axios";
import SolverResults from "./SolverResults";
import { useSearchParams } from "react-router-dom";

const GET_MOL_IMAGE_ENDPOINT = "/solver/get-mol-image";
const RUN_SOLVER_ENDPOINT = "/solver/run";

const Solver = () => {
  // before execution
  const [searchParams, setSearchParams] = useSearchParams({
    startingSmiles: "",
    targetSmiles: "",
  });

  const startingSmiles = searchParams.get("startingSmiles");
  const setStartingSmiles = useCallback(
    (smiles) => {
      setSearchParams((prev) => {
        prev.set("startingSmiles", smiles);
        return prev;
      });
    },
    [setSearchParams]
  );
  const [startingEncoding, setStartingEncoding] = useState(null);

  const targetSmiles = searchParams.get("targetSmiles");
  const setTargetSmiles = useCallback(
    (smiles) => {
      setSearchParams((prev) => {
        prev.set("targetSmiles", smiles);
        return prev;
      });
    },
    [setSearchParams]
  );
  const [targetEncoding, setTargetEncoding] = useState(null);

  // after execution
  const [solverResults, setSolverResults] = useState(null);

  // Effect to update state when searchParams change
  // Handles the display of the starting and target molecules,
  // both from the url and from the chemdraw in the PreExecution component
  useEffect(() => {
    const handleSetSolverMolecule = async (
      isStartingMolecule,
      preSolverSmiles
    ) => {
      const solverInfo = {
        true: {
          smiles: startingSmiles,
          setSmiles: setStartingSmiles,
          setEncoding: setStartingEncoding,
        },
        false: {
          smiles: targetSmiles,
          setSmiles: setTargetSmiles,
          setEncoding: setTargetEncoding,
        },
      };

      const { smiles, setSmiles, setEncoding } = solverInfo[isStartingMolecule];

      await axios
        .get(GET_MOL_IMAGE_ENDPOINT, {
          params: { smiles: smiles || preSolverSmiles },
        })
        .then((res) => {
          if (!smiles) {
            setSmiles(preSolverSmiles);
          }
          setEncoding(res.data);
        })
        .catch((error) => {
          alert(error.response.data.detail);
          setSmiles("");
        });
    };

    setSolverResults(null);

    if (startingSmiles && startingSmiles === targetSmiles) {
      alert("Starting and target molecules cannot be the same!");
      setTargetSmiles("");
      return;
    }

    if (startingSmiles) {
      handleSetSolverMolecule(true);
    } else {
      setStartingEncoding(null);
    }

    if (targetSmiles) {
      handleSetSolverMolecule(false);
    } else {
      setTargetEncoding(null);
    }
  }, [
    startingSmiles,
    targetSmiles,
    setStartingSmiles,
    setTargetSmiles,
    setStartingEncoding,
    setTargetEncoding,
  ]);

  // runs the solver and updates the solverResults state
  const handleRunSolver = async () => {
    await axios
      .get(RUN_SOLVER_ENDPOINT, {
        params: {
          start_smiles: startingSmiles,
          target_smiles: targetSmiles,
        },
      })
      .then((res) => {
        setSolverResults(res.data);
      });
  };

  if (solverResults) {
    return (
      <SolverResults
        startingSmiles={startingSmiles}
        targetSmiles={targetSmiles}
        solverResults={solverResults}
      />
    );
  } else {
    return (
      <PreExecution
        startingSmiles={startingSmiles}
        setStartingSmiles={setStartingSmiles}
        startingEncoding={startingEncoding}
        targetSmiles={targetSmiles}
        setTargetSmiles={setTargetSmiles}
        targetEncoding={targetEncoding}
        handleRunSolver={handleRunSolver}
      />
    );
  }
};

export default Solver;
