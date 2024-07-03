import axios from "axios";
import { useEffect, useState } from "react";
import { ALL_SUBSTRUCTURES_ENDPOINT } from "../../endpoints";

const ClassifierDocs = () => {
  const [allSubstructures, setAllSubstructures] = useState([]);
  useEffect(() => {
    const getAllSubstructures = async () => {
      await axios.get(ALL_SUBSTRUCTURES_ENDPOINT).then((res) => {
        setAllSubstructures(res.data);
      });
    };
    getAllSubstructures();
  }, []);

  return (
    <>
      <p>
        The classifier is a tool for identifying the substructures and
        functional groups within molecules. Start by drawing a molecule in the
        drawing tool on the left, and the substructures will be displayed on the
        right. The list of substructures will update in real time as the
        molecule is changed.
      </p>
      <img
        className="doc-image"
        src={process.env.PUBLIC_URL + "/docs/classifier-doc-main-page.png"}
        alt="Classifier main page doc"
      />
      <p>
        You can also open the current molecule in the Playground or the Solver
        by clicking on the options button.
      </p>
      <p>Here's the list of all available substructures:</p>
      {allSubstructures.length ? (
        <ul>
          {allSubstructures.map((substructure) => (
            <li key={substructure}>{substructure}</li>
          ))}
        </ul>
      ) : (
        <p>Loading...</p>
      )}
    </>
  );
};

export default ClassifierDocs;
