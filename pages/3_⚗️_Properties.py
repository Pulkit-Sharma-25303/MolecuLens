"""
Properties Calculator Page
Calculate and display molecular properties and descriptors.
"""
import streamlit as st
import pandas as pd
import sys
import os

# Add modules directory to sys.path for import
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from modules.molecular_utils import smiles_to_mol
from modules.descriptor_calc import get_mordred_df
from modules.screening_utils import lipinski_flags

st.set_page_config(page_title="Molecular Properties", page_icon="⚗️", layout="wide")

st.title("⚗️ Molecular Properties Calculator")
st.markdown("Calculate 2D/3D molecular descriptors, Lipinski filters, and physicochemical properties.")

# Sidebar input options
st.sidebar.header("Input Options")
input_source = st.sidebar.selectbox(
    "Input Source",
    ["Enter SMILES", "Upload File", "Use Example Molecule"]
)

smiles_list = []
data_ready = False

if input_source == "Enter SMILES":
    smiles_input = st.sidebar.text_area(
        "Paste SMILES (one per line)", height=80, value="CCO\nCC(N)C(=O)O"
    )
    if st.sidebar.button("Analyze"):
        smiles_list = [s.strip() for s in smiles_input.splitlines() if s.strip()]
        # Filter only valid SMILES
        smiles_list = [s for s in smiles_list if smiles_to_mol(s) is not None]
        data_ready = len(smiles_list) > 0

elif input_source == "Upload File":
    uploaded_file = st.sidebar.file_uploader("Upload CSV", type=["csv"])
    if uploaded_file and st.sidebar.button("Analyze"):
        df = pd.read_csv(uploaded_file)
        if 'smiles' in df.columns:
            smiles_list = [str(s).strip() for s in df['smiles'].tolist() if str(s).strip()]
            smiles_list = [s for s in smiles_list if smiles_to_mol(s) is not None]
            data_ready = len(smiles_list) > 0
        else:
            st.error("CSV file must contain a 'smiles' column.")

elif input_source == "Use Example Molecule":
    smiles_list = ["CCO", "CC(N)C(=O)O", "c1ccccc1"]
    data_ready = True

if data_ready and len(smiles_list) > 0:
    # Descriptor calculation
    desc_df = get_mordred_df(smiles_list)
    desc_df['Lipinski_Pass'] = [lipinski_flags(s) for s in smiles_list]
    desc_df['SMILES'] = smiles_list

    st.subheader("Properties Table")
    st.dataframe(desc_df, use_container_width=True)

    # Quick stats
    st.markdown("### Quick Stats")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Num Compounds", len(smiles_list))
    with col2:
        st.metric("Lipinski Pass", int(desc_df['Lipinski_Pass'].sum()))
    with col3:
        st.metric("Descriptors", len(desc_df.columns) - 2)  # Exclude Lipinski_Pass and SMILES

    # Download
    csv = desc_df.to_csv(index=False)
    st.download_button("Download as CSV", csv, "properties_results.csv", "text/csv")
else:
    st.info("Please provide SMILES or upload a supported file (CSV with a 'smiles' column) and click Analyze to calculate molecular properties.")

# Help section always visible
with st.expander("ℹ️ Properties Help"):
    st.markdown("""
    - **Molecular Descriptors**: Numerical features characterizing molecules.
    - **Lipinski's Rule of Five**: Checks drug-likeness properties.
    - **How to use**: Paste SMILES (one per line) or upload a CSV with a 'smiles' column; then click Analyze.
    """)
