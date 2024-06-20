import { Route, Routes } from "react-router-dom";
import Navbar from "./components/NavBar";
import Playground from "./pages/Playground";
import Solver from "./pages/Solver";
import Classifier from "./pages/Classifier";

const App = () => {
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
