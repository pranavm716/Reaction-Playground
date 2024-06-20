import { Route, Routes } from "react-router-dom";
import ChemDraw from "./ChemDraw";
import Navbar from "./NavBar";
import Playground from "./Playground";
import Solver from "./Solver";
import Classifier from "./Classifier";

export default function App() {
  return (
    <>
      <Navbar />
      <Routes>
        <Route path="/" element={<Playground />} />
        <Route path="/solver" element={<Solver />} />
        <Route path="/classifier" element={<Classifier />} />
      </Routes>
    </>
  );
}
