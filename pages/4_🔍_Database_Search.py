"""
Database Search Page
Search ChEMBL and PubChem databases
"""
import streamlit as st
import pandas as pd
import sys
import os

# Add modules to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from modules.database_client import search_chembl_by_target, search_pubchem_by_name, similarity_search
from modules.molecular_utils import smiles_to_mol

st.set_page_config(
    page_title="Database Search",
    page_icon="🔍",
    layout="wide"
)

st.title("🔍 Chemical Database Search")
st.markdown("Search ChEMBL bioactivity data and PubChem chemical information")

# Sidebar options
st.sidebar.header("Search Options")
database = st.sidebar.selectbox(
    "Database",
    ["ChEMBL", "PubChem", "Both"]
)

search_type = st.sidebar.selectbox(
    "Search Type",
    ["Target Search", "Compound Name", "SMILES Similarity", "Substructure"]
)

# Main search interface
st.subheader("Search Parameters")

col1, col2 = st.columns([2, 1])

with col1:
    if search_type == "Target Search":
        target_query = st.text_input(
            "Enter target name or UniProt ID:",
            placeholder="e.g., EGFR, P00533, Epidermal Growth Factor Receptor"
        )

        # Additional filters
        st.subheader("Filters")
        col3, col4 = st.columns(2)

        with col3:
            activity_type = st.selectbox(
                "Activity Type",
                ["Any", "IC50", "Ki", "EC50", "Kd"]
            )
            min_activity = st.number_input("Min Activity (nM)", min_value=0.0, value=0.0)

        with col4:
            max_activity = st.number_input("Max Activity (nM)", min_value=0.0, value=10000.0)
            confidence_score = st.slider("Confidence Score", 0, 9, 6)

    elif search_type == "Compound Name":
        compound_name = st.text_input(
            "Enter compound name:",
            placeholder="e.g., aspirin, caffeine, imatinib"
        )

    elif search_type == "SMILES Similarity":
        query_smiles = st.text_input(
            "Enter query SMILES:",
            placeholder="e.g., CCO"
        )
        similarity_threshold = st.slider("Similarity Threshold", 0.0, 1.0, 0.7, 0.05)

    elif search_type == "Substructure":
        substructure_smarts = st.text_input(
            "Enter SMARTS pattern:",
            placeholder="e.g., c1ccccc1 (benzene ring)"
        )

with col2:
    st.subheader("Quick Examples")

    if search_type == "Target Search":
        if st.button("🎯 EGFR"):
            target_query = "EGFR"
        if st.button("🎯 HIV Protease"):
            target_query = "HIV Protease"
        if st.button("🎯 Cyclooxygenase"):
            target_query = "COX"

    elif search_type == "Compound Name":
        if st.button("💊 Aspirin"):
            compound_name = "aspirin"
        if st.button("💊 Caffeine"):
            compound_name = "caffeine"
        if st.button("💊 Imatinib"):
            compound_name = "imatinib"

# Search button
search_button = st.button("🔍 Search Database", type="primary", use_container_width=True)

# Results section
if search_button:
    with st.spinner("Searching database..."):
        try:
            if search_type == "Target Search" and 'target_query' in locals():
                results = search_chembl_by_target(
                    target_query,
                    activity_type=activity_type if activity_type != "Any" else None,
                    min_activity=min_activity,
                    max_activity=max_activity,
                    confidence_score=confidence_score
                )

            elif search_type == "Compound Name" and 'compound_name' in locals():
                results = search_pubchem_by_name(compound_name)

            elif search_type == "SMILES Similarity" and 'query_smiles' in locals():
                results = similarity_search(query_smiles, similarity_threshold)

            else:
                results = pd.DataFrame()  # Empty results for demo

            # Display results
            if not results.empty:
                st.success(f"Found {len(results)} results")

                # Summary statistics
                st.subheader("Search Summary")
                col5, col6, col7, col8 = st.columns(4)

                with col5:
                    st.metric("Total Results", len(results))
                with col6:
                    if 'activity_value' in results.columns:
                        avg_activity = results['activity_value'].mean()
                        st.metric("Avg Activity", f"{avg_activity:.1f} nM")
                with col7:
                    if 'molecular_weight' in results.columns:
                        avg_mw = results['molecular_weight'].mean()
                        st.metric("Avg MW", f"{avg_mw:.1f} Da")
                with col8:
                    unique_compounds = results['canonical_smiles'].nunique() if 'canonical_smiles' in results.columns else len(results)
                    st.metric("Unique Compounds", unique_compounds)

                # Results table
                st.subheader("Results Table")

                # Add selection for columns to display
                all_cols = results.columns.tolist()
                default_cols = [col for col in ['compound_name', 'canonical_smiles', 'activity_value', 'activity_type', 'target_name'] if col in all_cols]

                selected_cols = st.multiselect(
                    "Select columns to display:",
                    all_cols,
                    default=default_cols[:5] if len(default_cols) >= 5 else default_cols
                )

                if selected_cols:
                    display_df = results[selected_cols].copy()
                    st.dataframe(display_df, use_container_width=True)

                # Export options
                st.subheader("Export Results")
                col9, col10, col11 = st.columns(3)

                with col9:
                    csv = results.to_csv(index=False)
                    st.download_button(
                        "📥 Download CSV",
                        csv,
                        "search_results.csv",
                        "text/csv",
                        use_container_width=True
                    )

                with col10:
                    if 'canonical_smiles' in results.columns:
                        sdf_data = "Mock SDF data"  # Would generate real SDF
                        st.download_button(
                            "📥 Download SDF",
                            sdf_data,
                            "compounds.sdf",
                            "chemical/x-mdl-sdfile",
                            use_container_width=True
                        )

                with col11:
                    st.button(
                        "🧪 Send to Editor",
                        use_container_width=True,
                        help="Send selected compounds to Molecular Editor"
                    )

            else:
                st.warning("No results found. Try adjusting your search parameters.")

        except Exception as e:
            st.error(f"Search error: {str(e)}")

# Database information
st.sidebar.subheader("Database Info")
st.sidebar.markdown("""
**ChEMBL:**
- 2.3M+ compounds
- 19M+ bioactivities
- 1.4M+ assays
- Manual curation

**PubChem:**
- 100M+ compounds
- Chemical properties
- Bioassay data
- Literature links
""")

# Help section
with st.expander("ℹ️ Search Help"):
    st.markdown("""
    ### Search Tips:

    **Target Search:**
    - Use gene names (EGFR, TP53)
    - UniProt IDs (P00533)
    - Protein names (Epidermal Growth Factor Receptor)

    **Compound Names:**
    - Common names (aspirin, caffeine)
    - IUPAC names
    - Trade names

    **SMILES Similarity:**
    - Use canonical SMILES
    - Tanimoto similarity threshold
    - Higher threshold = more similar

    **Substructure:**
    - SMARTS patterns
    - c1ccccc1 = benzene ring
    - N = any nitrogen
    """)

# Sample data for demonstration
if not search_button:
    st.subheader("Sample Results")
    st.info("Click search to see real results. Here's a sample of what you might find:")

    sample_data = {
        'Compound Name': ['Imatinib', 'Dasatinib', 'Nilotinib'],
        'SMILES': ['CCc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc4nccc(-c5cccnc5)n4',
                   'CNc1nc(Nc2ncc(s3)cc23)nc(N)c1C#N',
                   'CCc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc4nccc(-c5cccnc5)n4'],
        'Target': ['BCR-ABL', 'BCR-ABL', 'BCR-ABL'],
        'IC50 (nM)': [25, 0.8, 20],
        'MW (Da)': [493.6, 488.0, 529.5]
    }

    sample_df = pd.DataFrame(sample_data)
    st.dataframe(sample_df, use_container_width=True)
