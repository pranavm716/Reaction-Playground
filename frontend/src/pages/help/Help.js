import { useAppDoc } from "../../hooks";

const Help = ({ colorsMap }) => {
  const appDoc = useAppDoc(colorsMap);

  return (
    <div style={{ display: "flex", flexDirection: "column", gap: "15px" }}>
      <h2>Welcome to Reaction Playground!</h2>
      <p>
        This is a website intended to help beginner level organic chemistry
        students understand common organic reactions and functional groups.
        Reaction Playground consists of 3 main apps: The{" "}
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
    </div>
  );
};

export default Help;
