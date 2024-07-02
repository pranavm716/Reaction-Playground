import { useState } from "react";
import { NavLink, useNavigate } from "react-router-dom";

const appLinks = {
  "/": "Playground",
  "/solver": "Solver",
  "/classifier": "Classifier",
};

const docLinks = {
  "/help": "Help",
  "/about": "About",
};

const Navbar = ({ colorsMap }) => {
  const [hoverLink, setHoverLink] = useState(null);
  const navigate = useNavigate();

  const handleNav = (link, event) => {
    event.preventDefault();
    navigate(link);
  };

  const getNavLink = (link, text) => {
    return (
      <NavLink
        to={link}
        key={link}
        onMouseEnter={() => setHoverLink(link)}
        onMouseLeave={() => setHoverLink(null)}
        style={{ color: hoverLink === link ? colorsMap[link] : "#333" }}
        onClick={(event) => handleNav(link, event)}
      >
        {text}
      </NavLink>
    );
  };

  return (
    <nav className="navbar">
      <h1>Reaction Playground</h1>
      <div className="links">
        {Object.keys(appLinks).map((link) => {
          return getNavLink(link, appLinks[link]);
        })}

        <span style={{ padding: "0 10px", color: "#333" }}>|</span>

        {Object.keys(docLinks).map((link) => {
          return getNavLink(link, docLinks[link]);
        })}
      </div>
    </nav>
  );
};

export default Navbar;
