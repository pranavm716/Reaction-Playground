import axios from "axios";
import { useEffect, useState } from "react";
import { DocImageWithCaption } from "../../components/SmallUIComponents";
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
        functional groups within molecules.
      </p>
      <p>
        Start by drawing a molecule in the drawing tool on the left, and the
        substructures will be displayed on the right. The list of substructures
        will update in real time as the molecule is changed.
      </p>
      <DocImageWithCaption
        src="classifier-main-page.png"
        alt="Classifier main page doc"
        caption="Classifier page, with drawing tool on the left and substructures displayed on the right"
      />
      <p>
        You can also open the current molecule in the Playground or the Solver
        by clicking on the "Options" button:
      </p>
      <DocImageWithCaption
        src="classifier-options-button.png"
        alt="Classifier options button"
        caption="Options button in the Classifier tool"
      />
      <p>Here's the list of all available substructures:</p>
      {allSubstructures.length ? (
        <ul className="info-list">
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
