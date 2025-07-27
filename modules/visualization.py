"""
Molecular visualization utilities
"""

import logging
from typing import Optional, Union, Dict, Any

try:
    import stmol
    import py3Dmol
    STMOL_AVAILABLE = True
except ImportError:
    STMOL_AVAILABLE = False

try:
    from rdkit import Chem
    from rdkit.Chem import rdDepictor, Draw
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

def create_3d_view(mol, 
                  style: str = 'stick',
                  show_hydrogens: bool = True,
                  width: int = 400,
                  height: int = 400,
                  background: str = 'white') -> str:
    """
    Create 3D molecular visualization

    Args:
        mol: RDKit molecule object
        style: Visualization style ('stick', 'sphere', 'cartoon', 'surface')
        show_hydrogens: Whether to show hydrogen atoms
        width: Viewer width in pixels
        height: Viewer height in pixels
        background: Background color

    Returns:
        HTML string for 3D viewer
    """
    if not STMOL_AVAILABLE or not RDKIT_AVAILABLE:
        logging.warning("Required libraries not available for 3D visualization")
        return ""

    if mol is None:
        return ""

    try:
        # Add hydrogens if requested
        if show_hydrogens:
            mol = Chem.AddHs(mol)

        # Generate 3D coordinates
        from rdkit.Chem import AllChem
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol)

        # Convert to MOL block
        mol_block = Chem.MolToMolBlock(mol)

        # Create py3Dmol viewer
        viewer = py3Dmol.view(width=width, height=height)
        viewer.addModel(mol_block, 'mol')

        # Apply style
        if style == 'stick':
            viewer.setStyle({'stick': {'radius': 0.1}})
        elif style == 'sphere':
            viewer.setStyle({'sphere': {'radius': 0.3}})
        elif style == 'cartoon':
            viewer.setStyle({'cartoon': {}})
        elif style == 'surface':
            viewer.addSurface('VDW', {'opacity': 0.7})
            viewer.setStyle({'stick': {'radius': 0.1}})

        viewer.setBackgroundColor(background)
        viewer.zoomTo()

        return viewer._make_html()

    except Exception as e:
        logging.error(f"Error creating 3D visualization: {str(e)}")
        return ""

def create_2d_image(mol, 
                   size: tuple = (300, 300),
                   highlight_atoms: Optional[list] = None,
                   highlight_bonds: Optional[list] = None) -> bytes:
    """
    Create 2D molecular image

    Args:
        mol: RDKit molecule object
        size: Image size (width, height)
        highlight_atoms: List of atom indices to highlight
        highlight_bonds: List of bond indices to highlight

    Returns:
        Image bytes (PNG format)
    """
    if not RDKIT_AVAILABLE:
        return b""

    if mol is None:
        return b""

    try:
        # Generate 2D coordinates
        rdDepictor.Compute2DCoords(mol)

        # Create image
        img = Draw.MolToImage(
            mol, 
            size=size,
            highlightAtoms=highlight_atoms,
            highlightBonds=highlight_bonds
        )

        # Convert to bytes
        import io
        img_bytes = io.BytesIO()
        img.save(img_bytes, format='PNG')
        return img_bytes.getvalue()

    except Exception as e:
        logging.error(f"Error creating 2D image: {str(e)}")
        return b""

def create_grid_view(molecules: list,
                    mols_per_row: int = 4,
                    sub_img_size: tuple = (200, 200),
                    legends: Optional[list] = None) -> bytes:
    """
    Create grid view of multiple molecules

    Args:
        molecules: List of RDKit molecule objects
        mols_per_row: Number of molecules per row
        sub_img_size: Size of each molecule image
        legends: List of legends for each molecule

    Returns:
        Grid image bytes (PNG format)
    """
    if not RDKIT_AVAILABLE:
        return b""

    try:
        img = Draw.MolsToGridImage(
            molecules,
            molsPerRow=mols_per_row,
            subImgSize=sub_img_size,
            legends=legends
        )

        # Convert to bytes
        import io
        img_bytes = io.BytesIO()
        img.save(img_bytes, format='PNG')
        return img_bytes.getvalue()

    except Exception as e:
        logging.error(f"Error creating grid view: {str(e)}")
        return b""

def create_interactive_grid(smiles_list: list,
                          names: Optional[list] = None,
                          properties: Optional[Dict[str, list]] = None) -> str:
    """
    Create interactive molecular grid using mols2grid

    Args:
        smiles_list: List of SMILES strings
        names: List of molecule names
        properties: Dictionary of properties for each molecule

    Returns:
        HTML string for interactive grid
    """
    try:
        import mols2grid
        import pandas as pd

        # Prepare data
        data = {'SMILES': smiles_list}

        if names:
            data['Name'] = names
        else:
            data['Name'] = [f'Mol_{i:03d}' for i in range(len(smiles_list))]

        if properties:
            data.update(properties)

        df = pd.DataFrame(data)

        # Create grid
        grid = mols2grid.display(
            df,
            smiles_col='SMILES',
            name_col='Name',
            subset=['img'] + list(properties.keys()) if properties else ['img'],
            size=(150, 150)
        )

        return grid._repr_html_()

    except ImportError:
        logging.warning("mols2grid not available")
        return ""
    except Exception as e:
        logging.error(f"Error creating interactive grid: {str(e)}")
        return ""

def create_protein_view(pdb_string: str,
                       style: str = 'cartoon',
                       color_scheme: str = 'spectrum',
                       show_ligands: bool = True) -> str:
    """
    Create 3D protein visualization

    Args:
        pdb_string: PDB file content as string
        style: Visualization style
        color_scheme: Color scheme
        show_ligands: Whether to show ligands

    Returns:
        HTML string for 3D protein viewer
    """
    if not STMOL_AVAILABLE:
        return ""

    try:
        viewer = py3Dmol.view(width=600, height=400)
        viewer.addModel(pdb_string, 'pdb')

        # Protein visualization
        if style == 'cartoon':
            viewer.setStyle({'cartoon': {'color': color_scheme}})
        elif style == 'surface':
            viewer.addSurface('SAS', {'opacity': 0.7, 'color': color_scheme})
        elif style == 'ribbon':
            viewer.setStyle({'ribbon': {'color': color_scheme}})

        # Ligand visualization
        if show_ligands:
            viewer.setStyle({'hetflag': True}, {'stick': {'radius': 0.3}})

        viewer.zoomTo()
        return viewer._make_html()

    except Exception as e:
        logging.error(f"Error creating protein view: {str(e)}")
        return ""

def create_pharmacophore_view(mol,
                             pharmacophore_features: list) -> str:
    """
    Create pharmacophore visualization

    Args:
        mol: RDKit molecule object
        pharmacophore_features: List of pharmacophore features

    Returns:
        HTML string for pharmacophore viewer
    """
    # This would be more complex in a real implementation
    # For now, return basic 3D view
    return create_3d_view(mol, style='stick')

def create_reaction_view(reaction_smarts: str) -> str:
    """
    Create reaction visualization

    Args:
        reaction_smarts: SMARTS string representing the reaction

    Returns:
        HTML string for reaction viewer
    """
    if not RDKIT_AVAILABLE:
        return ""

    try:
        from rdkit.Chem import rdChemReactions

        rxn = rdChemReactions.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return ""

        # Create reaction image
        img = Draw.ReactionToImage(rxn, subImgSize=(300, 300))

        # Convert to base64 for HTML embedding
        import io
        import base64

        img_bytes = io.BytesIO()
        img.save(img_bytes, format='PNG')
        img_base64 = base64.b64encode(img_bytes.getvalue()).decode()

        html = f'<img src="data:image/png;base64,{img_base64}" alt="Reaction">'
        return html

    except Exception as e:
        logging.error(f"Error creating reaction view: {str(e)}")
        return ""

def create_similarity_map(reference_mol,
                         query_mol,
                         similarity_type: str = 'tanimoto') -> str:
    """
    Create molecular similarity visualization

    Args:
        reference_mol: Reference molecule
        query_mol: Query molecule
        similarity_type: Type of similarity calculation

    Returns:
        HTML string for similarity visualization
    """
    # This would require more complex implementation
    # For now, return side-by-side view
    try:
        # Create 2D images for both molecules
        ref_img = create_2d_image(reference_mol)
        query_img = create_2d_image(query_mol)

        # Create side-by-side HTML
        import base64

        ref_b64 = base64.b64encode(ref_img).decode()
        query_b64 = base64.b64encode(query_img).decode()

        html = f"""
        <div style="display: flex; justify-content: space-around;">
            <div>
                <h4>Reference</h4>
                <img src="data:image/png;base64,{ref_b64}" alt="Reference">
            </div>
            <div>
                <h4>Query</h4>
                <img src="data:image/png;base64,{query_b64}" alt="Query">
            </div>
        </div>
        """

        return html

    except Exception as e:
        logging.error(f"Error creating similarity map: {str(e)}")
        return ""

def export_visualization(viewer_html: str, 
                        filename: str,
                        format: str = 'html') -> bool:
    """
    Export visualization to file

    Args:
        viewer_html: HTML content of the viewer
        filename: Output filename
        format: Export format ('html', 'png', 'svg')

    Returns:
        Success status
    """
    try:
        if format == 'html':
            with open(filename, 'w') as f:
                f.write(viewer_html)
            return True
        else:
            logging.warning(f"Export format {format} not yet implemented")
            return False

    except Exception as e:
        logging.error(f"Error exporting visualization: {str(e)}")
        return False
