"""
MolecuLens - AI-Powered Drug Discovery Platform
Main Streamlit application entry point
"""
import streamlit as st
from pathlib import Path

# Configure page
st.set_page_config(
    page_title="MolecuLens",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        text-align: center;
        color: #2E86AB;
        padding: 2rem 0;
    }
    .feature-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 1rem;
        border-radius: 10px;
        color: white;
        margin: 1rem 0;
    }
    .sidebar .sidebar-content {
        background: linear-gradient(180deg, #2E86AB, #A23B72);
    }
</style>
""", unsafe_allow_html=True)

def main():
    # Sidebar navigation
    st.sidebar.image("assets/logo.png", width=200, use_column_width=True)
    st.sidebar.title("MolecuLens")
    st.sidebar.markdown("*AI-Powered Drug Discovery Platform*")

    # Main content
    st.markdown('<h1 class="main-header">🧬 Welcome to MolecuLens</h1>', unsafe_allow_html=True)

    st.markdown("""
    ## The Complete Drug Discovery Pipeline in Python

    MolecuLens is a comprehensive platform for computational drug discovery, built entirely with Python and Streamlit. 
    Explore molecular properties, search chemical databases, perform virtual screening, and simulate molecular dynamics - all through an intuitive web interface.
    """)

    # Feature overview
    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("""
        ### 🧪 Molecular Editor
        - Interactive structure drawing
        - SMILES input/output
        - 3D visualization
        - Real-time property calculation
        """)

    with col2:
        st.markdown("""
        ### ⚗️ Properties Calculator  
        - 1600+ molecular descriptors
        - Lipinski's Rule of Five
        - ADMET predictions
        - Batch processing
        """)

    with col3:
        st.markdown("""
        ### 🔍 Database Search
        - ChEMBL integration
        - PubChem queries
        - Bioactivity data
        - Similarity search
        """)

    col4, col5, col6 = st.columns(3)

    with col4:
        st.markdown("""
        ### 📊 Virtual Screening
        - Compound library filtering
        - Substructure search
        - Tanimoto similarity
        - Interactive results
        """)

    with col5:
        st.markdown("""
        ### 🎯 MD Simulation
        - OpenMM integration
        - Trajectory analysis
        - RMSD/RMSF plots
        - GPU acceleration
        """)

    with col6:
        st.markdown("""
        ### 🤖 AI Analytics
        - QSAR modeling
        - Machine learning
        - Predictive analysis
        - Pattern recognition
        """)

    # Getting started
    st.markdown("## 🚀 Getting Started")
    st.markdown("""
    1. **Navigate** using the sidebar to access different modules
    2. **Upload** your molecular data or draw structures
    3. **Analyze** properties and screen compounds
    4. **Simulate** molecular dynamics
    5. **Export** results for further analysis
    """)

    # Technical stack
    st.markdown("## 🛠️ Technology Stack")

    tech_col1, tech_col2 = st.columns(2)

    with tech_col1:
        st.markdown("""
        **Core Framework:**
        - Streamlit (Web UI)
        - Python 3.8+ (Backend)

        **Cheminformatics:**
        - RDKit (Molecular handling)
        - Mordred (Descriptors)
        - OpenEye (Advanced tools)
        """)

    with tech_col2:
        st.markdown("""
        **Visualization:**
        - Stmol (3D molecules)
        - Py3DMol (WebGL rendering)
        - Plotly (Interactive plots)

        **Simulation:**
        - OpenMM (Molecular dynamics)
        - MDAnalysis (Trajectory analysis)
        """)

if __name__ == "__main__":
    main()
