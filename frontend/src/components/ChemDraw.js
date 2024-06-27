import { Jsme } from "@loschmidt/jsme-react";

const ChemDraw = ({ smiles, setSmiles }) => {
  return (
    <div className="chem-draw-editor">
      <Jsme
        smiles={smiles}
        width={500}
        height={400}
        onChange={(smiles) => setSmiles(smiles)}
      />
    </div>
  );
};

export default ChemDraw;
