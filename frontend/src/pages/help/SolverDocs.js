import { DocImageWithCaption } from "../../components/SmallUIComponents";
import { useAllReactions } from "../../hooks";

const SolverDocs = () => {
  const { singleReactantReactions } = useAllReactions();

  return (
    <>
      <p>
        The solver is a tool for finding a synthetic reaction pathway
        between a starting and a target molecule.
      </p>

      <p>
        Start by drawing a starting molecule in the drawing tool on the left,
        then click the "Set as starting molecule" button to confirm your
        selection. Repeat the process for your target molecule. Here's what the
        page should look like after you've confirmed both molecules:
      </p>
      <DocImageWithCaption
        src="solver-pre-execution-page.png"
        alt="Solver pre-execution page"
        caption="Solver page, after the starting and target molecules are confirmed"
      />

      <p>
        If you need to clear the drawing canvas in between setting molecules,
        click the "Clear canvas" button located in the first row of buttons on
        the drawing tool:
      </p>
      <DocImageWithCaption
        src="solver-pre-execution-clear-canvas.png"
        alt="Solver pre-execution clear canvas"
        caption="Clear canvas button in the drawing tool"
      />

      <p>
        Then, click the "Find pathway" button to display the results. In some
        cases (like shown below), a reaction pathway is found. But in other
        cases, no such pathway may exist.
      </p>
      <DocImageWithCaption
        src="solver-results-page.png"
        alt="Solver results"
        caption="The Solver results page"
      />

      <p>
        You can hover over any reaction name to get a full description of what
        it does...
      </p>
      <DocImageWithCaption
        src="solver-reaction-description-tooltip.png"
        alt="Reaction description tooltip"
        caption="Tooltip displaying the description of a reaction"
      />

      <p>...and right click any molecule to open it in any of the apps.</p>
      <DocImageWithCaption
        src="solver-mol-right-click-menu.png"
        alt="Solver molecule right click menu"
        caption="Context menu for Solver molecules"
      />

      <p>
        <b style={{ color: "#b00020" }}>NOTE: </b>
        A limitation of the Solver is that it can only consider reactions that
        have a single reactant. This is because the Solver uses a brute force
        approach to find the reaction pathway and it has no way of knowing which
        additional reactants to choose that would push the products toward the
        target molecule. As a result, only the following single-reactant
        reactions are considered:
      </p>
      {singleReactantReactions.length ? (
        <ul className="info-list">
          {singleReactantReactions.map((reaction) => (
            <li key={reaction.name}>
              {reaction.name}
              <ul>
                <li>{reaction.description}</li>
              </ul>
            </li>
          ))}
        </ul>
      ) : (
        <p>Loading...</p>
      )}

      <p>
        <b style={{ color: "#b00020" }}>NOTE: </b>
        Another limitation is that reactions that produce multiple{" "}
        <em>products </em>
        may not yield the expected products when run on cyclic molecules. This
        is a work in progress.
      </p>
    </>
  );
};

export default SolverDocs;
