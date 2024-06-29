import React from "react";
import MolImage from "../../components/MolImage";
import {
  ArrowWithReactionInfo,
  plusIcon,
} from "../../components/SmallUIComponents";

const ProductPicker = ({
  products,
  setSmiles,
  reaction,
  molImage,
  missingReactantSmilesPicked,
  missingReactantEncodings,
}) => {
  let reactantMolRow;
  if (missingReactantSmilesPicked) {
    reactantMolRow = (
      <div className="mol-row">
        {molImage}
        {missingReactantSmilesPicked.map((smiles, index) => (
          <React.Fragment key={smiles}>
            {plusIcon}
            <MolImage
              smiles={smiles}
              encoding={missingReactantEncodings[index]}
            />
          </React.Fragment>
        ))}
      </div>
    );
  } else {
    reactantMolRow = molImage;
  }

  const pickProduct = (index) => {
    setSmiles(products[index].smiles);
  };

  return (
    <>
      {reactantMolRow}
      <ArrowWithReactionInfo
        reactionName={reaction.name}
        reactionDescriptionTooltip={reaction.description}
      />

      {/* list of clickable products */}
      <>
        <div className="mol-row">
          {products.map((product, index) => (
            <React.Fragment key={product.smiles}>
              <MolImage
                smiles={product.smiles}
                encoding={product.encoding}
                onClick={() => pickProduct(index)}
              />
              {index < products.length - 1 && plusIcon}
            </React.Fragment>
          ))}
        </div>
        <p className="multi-product-text">
          This produces {products.length} products. Click on the one you want to
          analyze next.
        </p>
      </>
    </>
  );
};

export default ProductPicker;
