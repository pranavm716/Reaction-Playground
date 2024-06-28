import { Tooltip } from "react-tooltip";
import React from "react";

const ReactionPicker = ({ reactions, setReactionPicked, molImage }) => {
  return (
    <>
      {molImage}
      {reactions.length ? (
        <>
          <p>Choose a reaction to run:</p>
          {reactions.map((reaction, index) => (
            <React.Fragment key={index}>
              <button
                className="reaction-button"
                onClick={() => setReactionPicked(reaction)}
              >
                <p>{reaction.name}</p>

                {/* info icon */}
                <svg
                  xmlns="http://www.w3.org/2000/svg"
                  height="24px"
                  viewBox="0 -960 960 960"
                  width="24px"
                  data-tooltip-id={`tt-${index}`}
                >
                  <path d="M450-290h60v-230h-60v230Zm30-298.46q13.73 0 23.02-9.29t9.29-23.02q0-13.73-9.29-23.02-9.29-9.28-23.02-9.28t-23.02 9.28q-9.29 9.29-9.29 23.02t9.29 23.02q9.29 9.29 23.02 9.29Zm.07 488.46q-78.84 0-148.21-29.92t-120.68-81.21q-51.31-51.29-81.25-120.63Q100-401.1 100-479.93q0-78.84 29.92-148.21t81.21-120.68q51.29-51.31 120.63-81.25Q401.1-860 479.93-860q78.84 0 148.21 29.92t120.68 81.21q51.31 51.29 81.25 120.63Q860-558.9 860-480.07q0 78.84-29.92 148.21t-81.21 120.68q-51.29 51.31-120.63 81.25Q558.9-100 480.07-100Zm-.07-60q134 0 227-93t93-227q0-134-93-227t-227-93q-134 0-227 93t-93 227q0 134 93 227t227 93Zm0-320Z" />
                </svg>
              </button>
              <Tooltip
                id={`tt-${index}`}
                place="bottom"
                offset={12}
                content={reaction.description}
                style={{ maxWidth: "350px", zIndex: 1000 }}
              />
            </React.Fragment>
          ))}
        </>
      ) : (
        <p>There are no compatible reactions.</p>
      )}
    </>
  );
};

export default ReactionPicker;
