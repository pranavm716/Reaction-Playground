import { Route, Routes, useLocation } from "react-router-dom";
import Navbar from "./components/NavBar";
import Playground from "./pages/playground/Playground";
import Solver from "./pages/solver/Solver";
import Classifier from "./pages/classifier/Classifier";
import { useEffect } from "react";

// TODO: Add help and about pages and color for classifier page
const primaryColorCSSTag = "--primary-color";

const useUpdatePrimaryColor = () => {
  const location = useLocation();

  useEffect(() => {
    const rootStyle = document.documentElement.style;
    switch (location.pathname) {
      case "/":
        rootStyle.setProperty(primaryColorCSSTag, "#f1356d");
        break;
      case "/solver":
        rootStyle.setProperty(primaryColorCSSTag, "#429987");
        break;
    }
  }, [location]);
}

const App = () => {
  useUpdatePrimaryColor();

  return (
    <>
      <Navbar />
      <div className="content">
        <Routes>
          <Route path="/" element={<Playground />} />
          <Route path="/solver" element={<Solver />} />
          <Route path="/classifier" element={<Classifier />} />
        </Routes>
      </div>
    </>
  );
}

export default App;
