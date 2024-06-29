import { useState, useEffect } from "react";
import { useSearchParams } from "react-router-dom";
import ChemDraw from "../../components/ChemDraw";
import { CautionText, MolImageMenu } from "../../components/SmallUIComponents";
import { useClassifierSubstructures } from "../../hooks";

export const CLASSIFIER_ENDPOINT = "/classifier/";

const Classifier = () => {
  const [searchParams, setSearchParams] = useSearchParams();
  const smiles = searchParams.get("smiles") || "";
  const { substructures, error } = useClassifierSubstructures(smiles);

  const updateSearchParams = (smiles) => {
    setSearchParams({ smiles });
  };

  const [chemDraw, setChemDraw] = useState(
    <ChemDraw setSmiles={updateSearchParams} />
  );

  // fill ChemDraw with smiles from URL on initial page load, if any
  useEffect(() => {
    if (smiles) {
      setChemDraw(<ChemDraw smiles={smiles} setSmiles={updateSearchParams} />);
    }

    // * This is needed due to a quirk in the JSME editor - we cannot set initial smiles and update the smiles from the same source (URL),
    // * so the dependencies array has to be empty
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  return (
    <>
      <p>
        This is the Classifier. Here, you can draw a molecule and find out which
        functional groups and substructures it contains.
      </p>
      <div className="two-panel-content">
        <div
          style={{
            display: "flex",
            flexDirection: "column",
            alignItems: "center",
          }}
        >
          {chemDraw}
          {smiles && <MolImageMenu smiles={smiles} />}
        </div>
        <div className="substructure-container">
          <p style={{ display: "flex", justifyContent: "center" }}>
            <b>Substructures</b>
          </p>
          <div className="substructures-display">
            {error ? (
              <CautionText text={error} />
            ) : substructures.length ? (
              substructures.map((substructure) => {
                return (
                  <button
                    disabled
                    className="substructure-chip"
                    key={substructure}
                  >
                    {substructure}
                  </button>
                );
              })
            ) : smiles ? (
              <p>No substructures detected.</p>
            ) : (
              <p>Draw a molecule to see its substructures!</p>
            )}
          </div>
        </div>
      </div>
    </>
  );
};

export default Classifier;
