import axios from "axios";
import React from "react";
import ReactDOM from "react-dom/client";
import Modal from "react-modal";
import { BrowserRouter } from "react-router-dom";
import App from "./App";
import "./index.css";

axios.defaults.withCredentials = true;
if (process.env.REACT_APP_IS_CLOUD_DEPLOYMENT === "true") {
  axios.defaults.headers.common["x-api-key"] = process.env.REACT_APP_API_KEY;
  axios.defaults.baseURL = process.env.REACT_APP_API_BASE_URL;
} else {
  axios.defaults.baseURL = "http://localhost:8000";
}

const root = ReactDOM.createRoot(document.getElementById("root"));
Modal.setAppElement("#root");

root.render(
  <React.StrictMode>
    <BrowserRouter>
      <App />
    </BrowserRouter>
  </React.StrictMode>,
);
