import { NavLink } from "react-router-dom";

const Navbar = () => {
    return (
        <nav className="navbar">
            <h1>Reaction Playground</h1>
            <div className="links">
                <NavLink to="/" >Playground</NavLink>
                <NavLink to="/solver">Solver</NavLink>
                <NavLink to="/classifier">Classifier</NavLink>
            </div>
        </nav>
    );
}

export default Navbar;