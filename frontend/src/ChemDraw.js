import { Jsme } from "@loschmidt/jsme-react";

const ChemDraw = () => {
    const logSmiles = (smiles) => {
        console.log(smiles)
    }
    return (
        <Jsme width={600} height={500} onChange={logSmiles}/>
    )
}
 
export default ChemDraw;