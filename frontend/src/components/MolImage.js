import { MenuItem, SubMenu, ControlledMenu } from '@szhsin/react-menu';
import '@szhsin/react-menu/dist/index.css';
import '@szhsin/react-menu/dist/transitions/slide.css';
import { useState } from 'react';

const openExternalIcon = <svg style={{
    marginLeft: '15px',
    position: 'relative',
}} xmlns="http://www.w3.org/2000/svg" height="20px" viewBox="0 -960 960 960" width="24px" fill="#000000">
    <path d="M212.31-140Q182-140 161-161q-21-21-21-51.31v-535.38Q140-778 161-799q21-21 51.31-21h252.3v60h-252.3q-4.62 0-8.46 3.85-3.85 3.84-3.85 8.46v535.38q0 4.62 3.85 8.46 3.84 3.85 8.46 3.85h535.38q4.62 0 8.46-3.85 3.85-3.84 3.85-8.46v-252.3h60v252.3Q820-182 799-161q-21 21-51.31 21H212.31Zm176.46-206.62-42.15-42.15L717.85-760H560v-60h260v260h-60v-157.85L388.77-346.62Z" />
</svg>

const copyIcon = <svg style={{
    marginLeft: '15px',
    position: 'relative',
}} xmlns="http://www.w3.org/2000/svg" height="20px" viewBox="0 -960 960 960" width="20px" fill="#000000">
    <path d="M362.31-260q-27.01 0-45.66-18.65Q298-297.3 298-324.31v-455.38q0-27.01 18.65-45.66Q335.3-844 362.31-844h359.38q27.01 0 45.66 18.65Q786-806.7 786-779.69v455.38q0 27.01-18.65 45.66Q748.7-260 721.69-260H362.31Zm0-52h359.38q4.62 0 8.46-3.85 3.85-3.84 3.85-8.46v-455.38q0-4.62-3.85-8.46-3.84-3.85-8.46-3.85H362.31q-4.62 0-8.46 3.85-3.85 3.84-3.85 8.46v455.38q0 4.62 3.85 8.46 3.84 3.85 8.46 3.85Zm-124 176q-27.01 0-45.66-18.65Q174-173.3 174-200.31v-507.38h52v507.38q0 4.62 3.85 8.46 3.84 3.85 8.46 3.85h411.38v52H238.31ZM350-312v-480 480Z" />
</svg>

const useContextMenu = (smiles) => {
    const [isOpen, setIsOpen] = useState(false);
    const [anchorPoint, setAnchorPoint] = useState({ x: 0, y: 0 });

    const contextMenu = <ControlledMenu
        anchorPoint={anchorPoint}
        state={isOpen ? 'open' : 'closed'}
        onClose={() => setIsOpen(false)}
        transition
    >
        <MenuItem onClick={() => navigator.clipboard.writeText(smiles)}>
            Copy SMILES
            {copyIcon}
        </MenuItem>
        <MenuItem>
            Open in Playground
            {openExternalIcon}
        </MenuItem>
        <SubMenu label='Open in Solver'>
            <MenuItem>
                As starting molecule
                {openExternalIcon}
            </MenuItem>
            <MenuItem>
                As target molecule
                {openExternalIcon}
            </MenuItem>
        </SubMenu>
        <MenuItem>
            Open in Classifier
            {openExternalIcon}
        </MenuItem>
    </ControlledMenu>;

    return { contextMenu, setIsOpen, setAnchorPoint };

};

const MolImage = ({ smiles, encoding, setSmiles, isHighlighted }) => {
    // Supports both clickable and non clickable mol images
    const className = setSmiles ? "clickable-mol-image" : "mol-image";
    const { contextMenu, setIsOpen, setAnchorPoint } = useContextMenu(smiles);

    return (
        <>
            <img
                onClick={setSmiles ? () => setSmiles(smiles) : undefined}
                className={className}
                width={200}
                height={200}
                src={`data:image/png;base64,${encoding}`}
                alt={smiles}
                style={isHighlighted ? { boxShadow: '0 0 8px rgba(0, 0, 0, 0.3)' } : undefined}
                onContextMenu={(e) => {
                    if (typeof document.hasFocus === 'function' && !document.hasFocus()) return;

                    e.preventDefault();
                    setAnchorPoint({ x: e.clientX, y: e.clientY });
                    setIsOpen(true);
                }}
            />
            {contextMenu}
        </>
    );

    // TODO: Make routes use url params and route context menu items to routes
}

export default MolImage;