import { ControlledMenu } from "@szhsin/react-menu";
import React, { useEffect, useState } from "react";
import { MolImageMenuItems, closeIcon } from "./components/SmallUIComponents";
import { useNavigate } from "react-router-dom";
import "@szhsin/react-menu/dist/index.css";
import "@szhsin/react-menu/dist/transitions/slide.css";
import axios from "axios";
import { CLASSIFIER_ENDPOINT } from "./endpoints";
import Modal from "react-modal";
import ChemDraw from "./components/ChemDraw";

export const useMolImageMenuReroute = (url, searchParams) => {
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

export const useClassifierSubstructures = (smiles) => {
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

export const useExtraReactantModal = (
  missingReactantSmiles,
  missingReactantPrompts,
  reactionName,
  setMissingReactantSmilesPicked
) => {
  const [modalIsOpen, setModalIsOpen] = useState(false);
  const [smiles, setSmiles] = useState(""); // Smiles of the current missing reactant picked

  // Modal styles and state management
  const customStyles = {
    content: {
      top: "50%",
      left: "50%",
      right: "auto",
      bottom: "auto",
      width: "60%",
      transform: "translate(-50%, -50%)",
    },
    closeIcon: {
      display: "flex",
      justifyContent: "flex-end",
    },
    header: {
      display: "flex",
      justifyContent: "space-between",
      alignItems: "center",
    },
    footer: {
      display: "flex",
      justifyContent: "center",
    },
  };

  const openModal = () => {
    setModalIsOpen(true);
  };

  const closeModal = () => {
    setModalIsOpen(false);
  };

  const handleUpdateMissingSmiles = (index) => {
    closeModal();
    missingReactantSmiles[index] = smiles;

    // If all missing reactants have been provided,
    // then update the state and let the parent's useEffect handle it from here
    if (missingReactantSmiles.every((smile) => smile !== null)) {
      setMissingReactantSmilesPicked(missingReactantSmiles);
    }
  };

  const modals = missingReactantPrompts.map((prompt, index) => (
    <Modal
      isOpen={modalIsOpen}
      onRequestClose={closeModal}
      contentLabel="Chem Draw Modal"
      style={customStyles}
    >
      <div style={customStyles.header}>
        {prompt} for {reactionName}
        <div style={customStyles.closeIcon}>
          {React.cloneElement(closeIcon, { onClick: closeModal })}
        </div>
      </div>

      <ChemDraw setSmiles={setSmiles} />

      <div style={customStyles.footer}>
        <button
          onClick={() => {
            handleUpdateMissingSmiles(index);
          }}
          className="primary-colored-button"
        >
          Save
        </button>
      </div>
    </Modal>
  ));

  return { openModal, modals };
};
