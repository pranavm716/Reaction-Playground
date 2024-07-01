import React from "react";
import MolImage from "../../components/MolImage";
import {
  ArrowWithReactionInfo,
  plusIcon,
} from "../../components/SmallUIComponents";

const molImageSideLength = 150;

const HistoryViewer = ({ history }) => {
  console.log(history);
  const prevHistory = history.slice(0, history.length - 1);
  if (!prevHistory.length) {
    return (
      <div className="history-viewer">
        <p>There is no history.</p>
      </div>
    );
  }
  return (
    <div className="history-viewer">
      {prevHistory.toReversed().map((step, index) => (
        <React.Fragment key={index}>
          <div className="mol-row">
            {step.productMetadata.map((product, index) => (
              <React.Fragment key={product.smiles}>
                <MolImage
                  smiles={product.smiles}
                  encoding={product.encoding}
                  sideLength={molImageSideLength}
                  isHighlighted={
                    step.productMetadata.length > 1 &&
                    index === step.productPicked
                  }
                />
                {index < step.productMetadata.length - 1 && plusIcon}
              </React.Fragment>
            ))}
          </div>
          <ArrowWithReactionInfo
            reactionName={step.reactionPicked.name}
            missingReactantMetadata={step.missingReactantMetadata}
            reactionDescriptionTooltip={step.reactionPicked.description}
            isDownArrow={false}
          />
        </React.Fragment>
      ))}

      {/* show the initial molecule  */}
      <MolImage
        smiles={prevHistory[0].curMolMetadata.smiles}
        encoding={prevHistory[0].curMolMetadata.encoding}
        sideLength={molImageSideLength}
      />
    </div>
  );
};

export default HistoryViewer;
