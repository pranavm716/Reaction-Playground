import { DocImageWithCaption } from "../../components/SmallUIComponents";
import { useAllReactions } from "../../hooks";

const PlaygroundDocs = () => {
  const { singleReactantReactions, multipleReactantReactions } =
    useAllReactions();

  return (
    <>
      <p>
        The playground is a sandbox for drawing and applying common organic
        reactions to molecules.
      </p>

      <p>
        Start by drawing a molecule in the drawing tool on the left, then click
        the "Run reactions" button.
      </p>
      <DocImageWithCaption
        src="playground-pre-loop-page.png"
        alt="Playground pre loop page"
        caption="Initial Playground page, after a molecule is drawn"
      />

      <p>
        A list of reactions compatible with this molecule will be shown on the
        left. Click the reaction you want to run on this molecule:
      </p>
      <DocImageWithCaption
        src="playground-pick-reaction.png"
        alt="Playground pick reaction"
        caption="Picking a reaction"
      />

      <p>
        Let's say we pick the hydrolysis reaction. There are 2 products for this
        reaction, as shown below. When there are multiple products, you will
        need to select the one to analyze next. If there is only 1 product, it
        will automatically get selected for you.
      </p>
      <DocImageWithCaption
        src="playground-pick-product.png"
        alt="Playground pick product"
        caption="Picking a product"
      />

      <p>
        Let's say we pick the second product, the amine. There is only 1
        available reaction for this molecule - Amide synthesis from acid
        chloride and amines. This reaction (among a few others; see the bottom
        of this page) requires additional reactants. In this case you will need
        to draw an acid chloride to proceed with the reaction:
      </p>
      <DocImageWithCaption
        src="playground-add-additional-reactants.png"
        alt="Playground add additional reactant"
        caption="Specifying additional reactants"
      />

      <p>
        Be sure to draw a valid molecule for this reaction, otherwise you will
        not be able to proceed. As you continue to apply reactions, a history of
        the molecules and reactions you have selected so far will be displayed
        on the right side:
      </p>
      <DocImageWithCaption
        src="playground-history-view.png"
        alt="Playground history view"
        caption="Playground history view"
      />

      <p>
        Clicking the "Previous molecule" button will allow you to step back
        through the history:
      </p>
      <DocImageWithCaption
        src="playground-previous-molecule.png"
        alt="Playground previous molecule"
        caption="Stepping back through history"
      />

      <p>
        Also, you can hover over any reaction name to get a full description of
        what it does...
      </p>
      <div
        style={{
          display: "flex",
          justifyContent: "center",
          alignItems: "center",
          gap: "50px",
        }}
      >
        <DocImageWithCaption
          src="playground-reaction-picker-description-tooltip.png"
          alt="Playground reaction picker tooltip"
          caption="Reaction description tooltip (when picking a reaction)"
        />
        <DocImageWithCaption
          src="playground-history-view-description-tooltip.png"
          alt="Playground history view tooltip"
          caption="Reaction description tooltip (in the history view)"
        />
      </div>

      <p>...and right click any molecule to open it in any of the apps.</p>
      <DocImageWithCaption
        src="playground-mol-right-click-menu.png"
        alt="Playgroud molcule right click menu"
        caption="Context menu for Playground molecules"
      />

      <p>Here's the list of all available reactions:</p>
      {singleReactantReactions.length ? (
        <>
          <u>Reactions with a single reactant</u>
          <ul className="info-list">
            {singleReactantReactions.map((reaction) => (
              <li>
                {reaction.name}
                <ul>
                  <li>{reaction.description}</li>
                </ul>
              </li>
            ))}
          </ul>

          <br />

          <u>Reactions with multiple reactants</u>
          <ul className="info-list">
            {multipleReactantReactions.map((reaction) => (
              <li>
                {reaction.name}
                <ul>
                  <li>{reaction.description}</li>
                </ul>
              </li>
            ))}
          </ul>
        </>
      ) : (
        <p>Loading...</p>
      )}

      <p>
        <b style={{ color: "#b00020" }}>NOTE: </b>A limitation of the Playground
        is that reactions that produce multiple <em>products </em>
        may not yield the expected products when run on cyclic molecules. This
        is a work in progress.
      </p>
    </>
  );
};

export default PlaygroundDocs;
