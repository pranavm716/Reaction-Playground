import { Jsme } from "@loschmidt/jsme-react";

const ChemDraw = ({setSmiles}) => {
    return (
        <div className="chem-draw-editor">
            <Jsme width={500} height={400} onChange={(smiles) => setSmiles(smiles)}/>
        </div>
    )
}
 
export default ChemDraw;