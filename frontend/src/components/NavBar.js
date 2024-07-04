import React, { useEffect, useState } from "react";
import { NavLink, useNavigate } from "react-router-dom";
import { Tooltip } from "react-tooltip";

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

  useEffect(() => {
    if (!localStorage.getItem("visited")) {
      localStorage.setItem("visited", "true");
    }
  }, []);

  const handleNav = (link, event) => {
    event.preventDefault();
    navigate(link);
  };

  const getNavLink = (link, text) => {
    let navLink = (
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

    if (link === "/help" && !localStorage.getItem("visited")) {
      navLink = (
        <>
          {React.cloneElement(navLink, {
            "data-tooltip-id": "initial-help-tooltip",
          })}
          <Tooltip
            id="initial-help-tooltip"
            place="bottom-end"
            offset={12}
            content="Welcome to Reaction Playground! If this is your first time here, check out the help page."
            defaultIsOpen={true}
            style={{ maxWidth: "350px", padding: "20px", zIndex: 1000 }}
            variant="info"
            clickable={true}
          />
        </>
      );
    }

    return navLink;
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
