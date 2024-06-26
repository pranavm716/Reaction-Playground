import { useState, useEffect, useCallback } from 'react';
import { useSearchParams } from 'react-router-dom';
import ChemDraw from '../../components/ChemDraw';
import axios from 'axios';
import { CautionText } from '../../components/SmallUIComponents';

const CLASSIFIER_ENDPOINT = '/classifier/';

const useSubstructures = (smiles) => {
    const [substructures, setSubstructures] = useState([]); // [list of str]
    const [error, setError] = useState(null);

    const retrieveSubstructures = useCallback(async () => {
        if (!smiles) {
            setSubstructures([]);
            setError(null);
            return;
        }
        await axios.get(CLASSIFIER_ENDPOINT, { params: { smiles } })
            .then(res => {
                setError(null);
                setSubstructures(res.data);
            })
            .catch(error => {
                setError(error.response.data.detail);
                setSubstructures([]);
            })
    }, [smiles])

    useEffect(() => { retrieveSubstructures() }, [retrieveSubstructures]);

    return { substructures, error };

};

const Classifier = () => {
    const [searchParams, setSearchParams] = useSearchParams({smiles: ''});
    const smiles = searchParams.get('smiles');
    const { substructures, error } = useSubstructures(smiles);

    const updateSearchParams = (smiles) => {
        setSearchParams({ smiles });
    }

    const [chemDraw, setChemDraw] = useState(<ChemDraw setSmiles={updateSearchParams} />);
    useEffect(() => {
        if (smiles) {
            setChemDraw(<ChemDraw smiles={smiles} setSmiles={updateSearchParams} />);
        }
    }, []);

    return (
        <>
            <p>
                This is the Classifier. Here, you can draw a molecule and find out which functional groups and substructures it contains.
            </p>
            <div className="two-panel-content">
                <div>
                    {chemDraw}
                </div>
                <div className='substructure-container'>
                    <p style={{ display: 'flex', justifyContent: 'center' }}>
                        <b>Substructures</b>
                    </p>
                    <div className='substructures-display'>
                        {
                            error ?
                                <CautionText text={error} />
                                : substructures.length ?
                                    substructures.map(substructure => {
                                        return (
                                            <button disabled className="substructure-chip" key={substructure}>
                                                {substructure}
                                            </button>
                                        )
                                    })
                                    :
                                    smiles ?
                                        <p>No substructures detected.</p>
                                        :
                                        <p>Draw a molecule to see its substructures!</p>
                        }
                    </div>
                </div>
            </div>
        </>
    );
}

export default Classifier;