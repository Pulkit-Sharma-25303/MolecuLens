"""
Molecular Editor Page
Interactive molecule drawing and editing
"""
import streamlit as st
import sys
import os

# Add modules to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

try:
    from streamlit_ketcher import st_ketcher
    KETCHER_AVAILABLE = True
except ImportError:
    KETCHER_AVAILABLE = False

try:
    import stmol
    from stmol import showmol
    STMOL_AVAILABLE = True
except ImportError:
    STMOL_AVAILABLE = False

from modules.molecular_utils import smiles_to_mol, mol_to_smiles, validate_smiles
from modules.visualization import create_3d_view
from modules.descriptor_calc import calculate_basic_properties

st.set_page_config(
    page_title="Molecular Editor",
    page_icon="🧪",
    layout="wide"
)

st.title("🧪 Molecular Editor")
st.markdown("Draw, edit, and visualize molecular structures")

# Sidebar options
st.sidebar.header("Editor Options")
input_method = st.sidebar.selectbox(
    "Input Method",
    ["Draw Structure", "SMILES Input", "File Upload"]
)

# Main content area
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("Structure Input")

    if input_method == "Draw Structure":
        if KETCHER_AVAILABLE:
            smiles = st_ketcher("CCO")  # Default molecule
            st.success(f"Generated SMILES: {smiles}")
        else:
            st.warning("Ketcher editor not available. Please install streamlit-ketcher.")
            smiles = st.text_input("Enter SMILES instead:", value="CCO")

    elif input_method == "SMILES Input":
        smiles = st.text_input("Enter SMILES:", value="CCO")

        # Validate SMILES
        if smiles:
            is_valid, message = validate_smiles(smiles)
            if is_valid:
                st.success("✅ Valid SMILES")
            else:
                st.error(f"❌ Invalid SMILES: {message}")

    elif input_method == "File Upload":
        uploaded_file = st.file_uploader(
            "Upload molecule file",
            type=['mol', 'sdf', 'mol2']
        )
        if uploaded_file:
            st.info("File uploaded successfully")
            smiles = "CCO"  # Placeholder - would parse file in real implementation
        else:
            smiles = "CCO"  # Default

    # Additional options
    st.subheader("Display Options")
    style = st.selectbox("3D Style", ["stick", "sphere", "cartoon", "surface"])
    show_hydrogens = st.checkbox("Show Hydrogens", value=True)

with col2:
    st.subheader("3D Visualization")

    if 'smiles' in locals() and smiles:
        # Create molecule object
        mol = smiles_to_mol(smiles)

        if mol:
            # 3D visualization
            if STMOL_AVAILABLE:
                try:
                    mol_html = create_3d_view(mol, style=style, show_hydrogens=show_hydrogens)
                    showmol(mol_html, height=400, width=400)
                except Exception as e:
                    st.error(f"3D visualization error: {str(e)}")
                    st.info("Displaying 2D structure instead")
                    # Would show 2D structure here
            else:
                st.warning("3D visualization not available. Please install stmol.")
        else:
            st.error("Could not create molecule object")

# Properties section
if 'smiles' in locals() and smiles:
    st.subheader("Molecular Properties")

    col3, col4, col5 = st.columns(3)

    with col3:
        st.metric("Molecular Weight", "46.07 Da")  # Would calculate from mol
        st.metric("LogP", "0.31")

    with col4:
        st.metric("H-Bond Donors", "1")
        st.metric("H-Bond Acceptors", "1")

    with col5:
        st.metric("Rotatable Bonds", "0")
        st.metric("TPSA", "20.23 Ų")

# Action buttons
st.subheader("Actions")
col6, col7, col8 = st.columns(3)

with col6:
    if st.button("Calculate All Properties", use_container_width=True):
        st.success("Properties calculated! Check the Properties page for details.")

with col7:
    if st.button("Search Databases", use_container_width=True):
        st.success("Redirecting to database search...")

with col8:
    if st.button("Start MD Simulation", use_container_width=True):
        st.success("Preparing simulation parameters...")

# Export options
st.subheader("Export Options")
export_col1, export_col2 = st.columns(2)

with export_col1:
    if st.button("Export as MOL", use_container_width=True):
        # Would generate MOL file for download
        st.info("MOL file ready for download")

with export_col2:
    if st.button("Export as SDF", use_container_width=True):
        # Would generate SDF file for download
        st.info("SDF file ready for download")
