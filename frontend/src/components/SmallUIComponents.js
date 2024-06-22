import React from "react";

export const cautionIcon = <svg style={{
    paddingRight: '5px',
    position: 'relative',
    top: '5px',
}} xmlns="http://www.w3.org/2000/svg" height="24px" viewBox="0 -960 960 960" width="24px" fill="#f00">
    <path d="M109.23-160 480-800l370.77 640H109.23ZM178-200h604L480-720 178-200Zm302-55.38q10.46 0 17.54-7.08 7.08-7.08 7.08-17.54 0-10.46-7.08-17.54-7.08-7.08-17.54-7.08-10.46 0-17.54 7.08-7.08 7.08-7.08 17.54 0 10.46 7.08 17.54 7.08 7.08 17.54 7.08Zm-20-89.24h40v-200h-40v200ZM480-460Z" />
</svg>;

export const closeIcon = <svg style={{
    cursor: 'pointer',
}} xmlns="http://www.w3.org/2000/svg" height="20px" viewBox="0 -960 960 960" width="20px" fill="#000000">
    <path d="m252.62-217.23-35.39-35.39L444.62-480 217.23-707.38l35.39-35.39L480-515.38l227.38-227.39 35.39 35.39L515.38-480l227.39 227.38-35.39 35.39L480-444.62 252.62-217.23Z" />
</svg>;

export const plusIcon = <svg style={{
    padding: "0 10px",
}} xmlns="http://www.w3.org/2000/svg" height="24px" viewBox="0 -960 960 960" width="24px" fill="#000000">
    <path d="M450-450H220v-60h230v-230h60v230h230v60H510v230h-60v-230Z" />
</svg>;

export const ArrowWithReactionName = ({ reactionName }) => {
    return (
        <div style={{
            position: 'relative',
            display: 'flex',
            justifyContent: 'center',
            alignItems: 'center',
        }}>
            <img height={130} src={process.env.PUBLIC_URL + "down-arrow.png"}></img>
            <span style={{ width: 200, position: 'absolute', left: 60 }}>{reactionName}</span>
        </div>
    );
};

