"""
MolecuLens Home Page
"""
import streamlit as st
from pathlib import Path

st.set_page_config(
    page_title="MolecuLens - Home",
    page_icon="🏠",
    layout="wide"
)

st.title("🏠 MolecuLens Home")
st.markdown("### Welcome to the AI-Powered Drug Discovery Platform")

# Quick navigation
st.markdown("## Quick Navigation")

col1, col2, col3 = st.columns(3)

with col1:
    if st.button("🧪 Molecular Editor", use_container_width=True):
        st.switch_page("pages/2_🧪_Molecular_Editor.py")
    if st.button("📊 Virtual Screening", use_container_width=True):
        st.switch_page("pages/5_📊_Virtual_Screening.py")

with col2:
    if st.button("⚗️ Properties Calculator", use_container_width=True):
        st.switch_page("pages/3_⚗️_Properties.py")
    if st.button("🎯 MD Simulation", use_container_width=True):
        st.switch_page("pages/6_🎯_MD_Simulation.py")

with col3:
    if st.button("🔍 Database Search", use_container_width=True):
        st.switch_page("pages/4_🔍_Database_Search.py")

# Recent activity
st.markdown("## Recent Activity")
st.info("Welcome! Start by drawing a molecule in the Molecular Editor or searching the database.")

# Tips
st.markdown("## 💡 Tips for Getting Started")
st.markdown("""
- **New to drug discovery?** Start with the Molecular Editor to explore basic concepts
- **Have specific compounds?** Use the Database Search to find bioactivity data
- **Working with libraries?** Try Virtual Screening to filter large datasets
- **Need simulations?** The MD Simulation module supports OpenMM workflows
""")
