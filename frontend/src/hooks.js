import { ControlledMenu } from "@szhsin/react-menu";
import { useEffect, useState } from "react";
import { MolImageMenuItems } from "./components/SmallUIComponents";
import { useNavigate } from "react-router-dom";
import "@szhsin/react-menu/dist/index.css";
import "@szhsin/react-menu/dist/transitions/slide.css";
import axios from "axios";
import { CLASSIFIER_ENDPOINT } from "./pages/classifier/Classifier";

export const useReroute = (url, searchParams) => {
  const navigate = useNavigate();

  return () => {
    navigate(`${url}?${new URLSearchParams(searchParams).toString()}`);
  };
};

export const useMolImageContextMenu = (smiles) => {
  const [isOpen, setIsOpen] = useState(false);
  const [anchorPoint, setAnchorPoint] = useState({ x: 0, y: 0 });

  const contextMenu = (
    <ControlledMenu
      anchorPoint={anchorPoint}
      state={isOpen ? "open" : "closed"}
      onClose={() => setIsOpen(false)}
      transition
    >
      <MolImageMenuItems smiles={smiles} />
    </ControlledMenu>
  );

  return { contextMenu, setIsOpen, setAnchorPoint };
};

export const useSubstructures = (smiles) => {
  const [substructures, setSubstructures] = useState([]); // [list of str]
  const [error, setError] = useState(null);

  useEffect(() => {
    const retrieveSubstructures = async () => {
      if (!smiles) {
        setSubstructures([]);
        setError(null);
        return;
      }
      await axios
        .get(CLASSIFIER_ENDPOINT, { params: { smiles } })
        .then((res) => {
          setError(null);
          setSubstructures(res.data);
        })
        .catch((error) => {
          setError(error.response.data.detail);
          setSubstructures([]);
        });
    };

    retrieveSubstructures();
  }, [smiles]);

  return { substructures, error };
};
