import { useAppDoc } from "../../hooks";

const Help = ({ colorsMap }) => {
  const appDoc = useAppDoc(colorsMap);

  return (
    <div style={{ display: "flex", flexDirection: "column", gap: "15px" }}>
      <h2>Welcome to Reaction Playground!</h2>
      <p>
        This is a tool that allows users to experiment with applying common
        organic reactions to molecules. It is intended to help beginner level
        organic chemistry students. Reaction Playground consists of 3 main apps:
        The{" "}
        <span>
          <b style={{ color: colorsMap["/"] }}>Playground</b>
        </span>
        , the{" "}
        <span>
          <b style={{ color: colorsMap["/solver"] }}>Solver</b>
        </span>
        , and the{" "}
        <span>
          <b style={{ color: colorsMap["/classifier"] }}>Classifier</b>
        </span>
        .
      </p>
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
// TODO: Add a tooltip directing user to the help page on first load
