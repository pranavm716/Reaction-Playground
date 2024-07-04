const PlaygroundDocs = () => {
  return (
    <>
      <p>
        The playground is a sandbox for drawing and applying common organic
        reactions to molecules.
      </p>

      <p>
        <b style={{ color: "#b00020" }}>NOTE: </b>
        A limitation of the Playground is that reactions that produce multiple{" "}
        <em>products </em>
        may not yield the expected products when run on cyclic molecules. This
        is a work in progress.
      </p>
    </>
  );
};

export default PlaygroundDocs;
