import { NavLink, useNavigate, useLocation } from "react-router-dom";
import { useState } from "react";

const links = { "/": "Playground", "/solver": "Solver", "/classifier": "Classifier" }

const Navbar = ({ colorsMap }) => {
    const [hoverLink, setHoverLink] = useState(null);
    const navigate = useNavigate();
    const location = useLocation();

    const handleNav = (link, event) => {
        event.preventDefault();
        navigate(link);
        if (location.pathname === link) {
            window.location.reload();
        }
    };

    return (
        <nav className="navbar">
            <h1>Reaction Playground</h1>
            <div className="links">
                {Object.keys(links).map((link) => {
                    return (
                        <NavLink
                            to={link}
                            key={link}
                            onMouseEnter={() => setHoverLink(link)}
                            onMouseLeave={() => setHoverLink(null)}
                            style={{ color: hoverLink === link ? colorsMap[link] : "#333" }}
                            onClick={(event) => handleNav(link, event)}
                        >
                            {links[link]}
                        </NavLink>
                    );
                })}
            </div>
        </nav>
    );
}

export default Navbar;