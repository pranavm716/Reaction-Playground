import { useAppDoc } from "../../hooks";
const appInfo = {
  "/": {
    name: "Playground",
    description: "A sandbox for drawing and applying reactions to molecules.",
  },
  "/solver": {
    name: "Solver",
    description:
      "A tool for calculating a synthetic reaction pathway between a starting and a target molecule.",
  },
  "/classifier": {
    name: "Classifier",
    description:
      "A tool for identifying the substructures and functional groups within molecules.",
  },
};

const Help = ({ colorsMap }) => {
  const appDoc = useAppDoc(colorsMap);

  return (
    <div style={{ display: "flex", flexDirection: "column", gap: "15px" }}>
      <h2>Welcome to Reaction Playground!</h2>
      <p>
        This is a tool that allows users to experiment with applying common
        organic reactions to molecules. It is intended to help beginner organic
        chemistry students. Reaction Playground consists of 3 main apps:
      </p>
      <ul>
        {Object.keys(appInfo).map((app) => {
          return (
            <li key={app}>
              <span>
                <b style={{ color: colorsMap[app] }}>{appInfo[app].name}</b>
              </span>
              : {appInfo[app].description}
            </li>
          );
        })}
      </ul>
      <p>Click through the tabs below to learn more about each app.</p>
      {appDoc}
      {/* <ul>
        <li>Intramolecular reactions do not function properly.</li>
        <li>
          Use of the browser back and forward buttons may not function properly
          when: 1) Switching between routes/apps (you may need to click multiple
          times) 2) Within the classifier app, since it has a live update with
          the ChemDraw tool.
        </li>
      </ul> */}
    </div>
  );
};

export default Help;
