import { useEffect } from "react";
import { Route, Routes, useLocation } from "react-router-dom";
import Navbar from "./components/NavBar";
import Classifier from "./pages/classifier/Classifier";
import About from "./pages/docs/About";
import Help from "./pages/docs/Help";
import Playground from "./pages/playground/Playground";
import Solver from "./pages/solver/Solver";

// TODO: Add help and about pages
const colorsMap = {
  "/": "#f1356d",
  "/solver": "#429987",
  "/classifier": "#1976D2",
  "/help": "#7E57C2",
  "/about": "#DA552F",
};

const useUpdatePrimaryColor = () => {
  const location = useLocation();

  useEffect(() => {
    const primaryColor = colorsMap[location.pathname];
    document.documentElement.style.setProperty("--primary-color", primaryColor);
  }, [location]);
};

const App = () => {
  useUpdatePrimaryColor();

  return (
    <>
      <Navbar colorsMap={colorsMap} />
      <div className="content">
        <Routes>
          <Route path="/" element={<Playground />} />
          <Route path="/solver" element={<Solver />} />
          <Route path="/classifier" element={<Classifier />} />
          <Route path="/help" element={<Help />} />
          <Route path="/about" element={<About />} />
        </Routes>
      </div>
    </>
  );
};

export default App;
