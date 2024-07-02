const references = [
  "https://www.rdkit.org/docs/index.html",
  "https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html",
  "https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html",
  "https://github.com/rdkit/rdkit-tutorials/tree/master/notebooks",
  "https://chem.libretexts.org/",
  "https://stackabuse.com/courses/graphs-in-python-theory-and-implementation/lessons/breadth-first-search-bfs-algorithm/",
  "https://www.youtube.com/watch?v=j942wKiXFu8&list=PL4cUxeGkcC9gZD-Tvwfod2gaISzfRiP9d"
  // TODO: link for deploying app on GH pages
];

const About = () => {
  return (
    <>
      <div className="about-section">
        <h2>Source code</h2>
        <p>
          See the source code on GitHub{" "}
          <a href="https://github.com/pranavm716/Reaction-Playground">here</a>.
        </p>
      </div>

      <div className="about-section">
        <h2>Tech stack</h2>
        <ul>
          <li>
            Backend: <a href="https://www.rdkit.org/">RDKit</a>, a
            cheminformatics library (Python)
          </li>
          <li>
            Frontend: <a href="https://react.dev/">React</a>, a frontend UI
            framework (JavaScript)
          </li>
          <li>
            API layer: <a href="https://fastapi.tiangolo.com/">FastAPI</a>, a
            web framework for building REST APIs (Python)
          </li>
        </ul>
      </div>

      <div className="about-section">
        <h2>References</h2>
        <ul>
          {references.map((ref, index) => (
            <li key={index}>
              <a href={ref}>{ref}</a>
            </li>
          ))}
          <li>My own organic chemistry notes :)</li>
        </ul>
      </div>
    </>
  );
};

export default About;
