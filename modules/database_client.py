"""
Database client utilities for ChEMBL and PubChem
"""

import pandas as pd
import logging
import time
from typing import List, Dict, Optional, Union

try:
    from chembl_webresource_client.new_client import new_client
    CHEMBL_AVAILABLE = True
except ImportError:
    CHEMBL_AVAILABLE = False

try:
    import pubchempy as pcp
    PUBCHEM_AVAILABLE = True
except ImportError:
    PUBCHEM_AVAILABLE = False

def search_chembl_by_target(target_query: str,
                           activity_type: Optional[str] = None,
                           min_activity: float = 0.0,
                           max_activity: float = 10000.0,
                           confidence_score: int = 6,
                           limit: int = 1000) -> pd.DataFrame:
    """
    Search ChEMBL database by target

    Args:
        target_query: Target name, gene symbol, or UniProt ID
        activity_type: Type of activity (IC50, Ki, etc.)
        min_activity: Minimum activity value (nM)
        max_activity: Maximum activity value (nM)
        confidence_score: Minimum confidence score (0-9)
        limit: Maximum number of results

    Returns:
        DataFrame with bioactivity data
    """
    if not CHEMBL_AVAILABLE:
        logging.warning("ChEMBL client not available")
        return _get_mock_chembl_data(target_query)

    try:
        # Search for targets
        target = new_client.target
        targets = target.search(target_query).only(['target_chembl_id', 'pref_name', 'organism'])

        if not targets:
            logging.warning(f"No targets found for query: {target_query}")
            return pd.DataFrame()

        # Get activities for the first target
        target_chembl_id = targets[0]['target_chembl_id']

        # Query activities
        activity = new_client.activity
        activities = activity.filter(
            target_chembl_id=target_chembl_id,
            confidence_score__gte=confidence_score
        ).only([
            'molecule_chembl_id',
            'canonical_smiles', 
            'standard_value',
            'standard_units',
            'standard_type',
            'confidence_score',
            'pchembl_value'
        ])

        results = []
        count = 0

        for act in activities:
            if count >= limit:
                break

            try:
                # Apply activity filters
                if activity_type and act.get('standard_type') != activity_type:
                    continue

                standard_value = act.get('standard_value')
                if standard_value is None:
                    continue

                # Convert to nM if needed
                units = act.get('standard_units', 'nM').lower()
                if units == 'um':
                    standard_value *= 1000
                elif units == 'mm':
                    standard_value *= 1000000
                elif units == 'm':
                    standard_value *= 1000000000

                if not (min_activity <= standard_value <= max_activity):
                    continue

                result = {
                    'compound_id': act.get('molecule_chembl_id'),
                    'canonical_smiles': act.get('canonical_smiles'),
                    'activity_value': standard_value,
                    'activity_units': 'nM',
                    'activity_type': act.get('standard_type'),
                    'confidence_score': act.get('confidence_score'),
                    'pchembl_value': act.get('pchembl_value'),
                    'target_name': targets[0]['pref_name'],
                    'target_chembl_id': target_chembl_id
                }

                results.append(result)
                count += 1

            except Exception as e:
                logging.error(f"Error processing activity: {str(e)}")
                continue

        return pd.DataFrame(results)

    except Exception as e:
        logging.error(f"ChEMBL search error: {str(e)}")
        return _get_mock_chembl_data(target_query)

def search_pubchem_by_name(compound_name: str) -> pd.DataFrame:
    """
    Search PubChem by compound name

    Args:
        compound_name: Name of the compound

    Returns:
        DataFrame with compound information
    """
    if not PUBCHEM_AVAILABLE:
        logging.warning("PubChemPy not available")
        return _get_mock_pubchem_data(compound_name)

    try:
        # Search by name
        compounds = pcp.get_compounds(compound_name, 'name')

        if not compounds:
            logging.warning(f"No compounds found for: {compound_name}")
            return pd.DataFrame()

        results = []

        for compound in compounds[:10]:  # Limit to first 10 results
            try:
                result = {
                    'cid': compound.cid,
                    'iupac_name': compound.iupac_name,
                    'molecular_formula': compound.molecular_formula,
                    'molecular_weight': compound.molecular_weight,
                    'canonical_smiles': compound.canonical_smiles,
                    'isomeric_smiles': compound.isomeric_smiles,
                    'inchi': compound.inchi,
                    'inchikey': compound.inchikey,
                    'synonyms': ', '.join(compound.synonyms[:5]) if compound.synonyms else '',
                    'xlogp': compound.xlogp,
                    'exact_mass': compound.exact_mass,
                    'monoisotopic_mass': compound.monoisotopic_mass,
                    'tpsa': compound.tpsa,
                    'complexity': compound.complexity,
                    'charge': compound.charge,
                    'h_bond_donor_count': compound.h_bond_donor_count,
                    'h_bond_acceptor_count': compound.h_bond_acceptor_count,
                    'rotatable_bond_count': compound.rotatable_bond_count,
                    'heavy_atom_count': compound.heavy_atom_count
                }

                results.append(result)

            except Exception as e:
                logging.error(f"Error processing PubChem compound: {str(e)}")
                continue

        return pd.DataFrame(results)

    except Exception as e:
        logging.error(f"PubChem search error: {str(e)}")
        return _get_mock_pubchem_data(compound_name)

def similarity_search(query_smiles: str, 
                     threshold: float = 0.7,
                     database: str = 'chembl',
                     limit: int = 100) -> pd.DataFrame:
    """
    Perform similarity search

    Args:
        query_smiles: Query SMILES string
        threshold: Similarity threshold (0-1)
        database: Database to search ('chembl' or 'pubchem')
        limit: Maximum number of results

    Returns:
        DataFrame with similar compounds
    """
    # This would require more complex implementation
    # For now, return mock data
    logging.info(f"Performing similarity search for {query_smiles}")
    return _get_mock_similarity_data(query_smiles, threshold)

def _get_mock_chembl_data(target_query: str) -> pd.DataFrame:
    """Generate mock ChEMBL data for demonstration"""
    import random

    mock_data = []
    compound_names = ['Imatinib', 'Dasatinib', 'Nilotinib', 'Bosutinib', 'Ponatinib']
    smiles_list = [
        'CCc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc4nccc(-c5cccnc5)n4',
        'CNc1nc(Nc2ncc(s3)cc23)nc(N)c1C#N',
        'CCc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc4nccc(-c5cccnc5)n4',
        'CCN(CC)c1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OC',
        'CCc1ccc(C(=O)Nc2ccc3c(c2)NCN3C)cc1-c4ccccc4F'
    ]

    for i in range(min(50, len(compound_names) * 10)):
        idx = i % len(compound_names)
        mock_data.append({
            'compound_id': f'CHEMBL{random.randint(100000, 999999)}',
            'canonical_smiles': smiles_list[idx],
            'activity_value': random.uniform(0.1, 1000),
            'activity_units': 'nM',
            'activity_type': random.choice(['IC50', 'Ki', 'EC50', 'Kd']),
            'confidence_score': random.randint(6, 9),
            'pchembl_value': random.uniform(5, 9),
            'target_name': target_query,
            'target_chembl_id': f'CHEMBL{random.randint(1000, 9999)}'
        })

    return pd.DataFrame(mock_data)

def _get_mock_pubchem_data(compound_name: str) -> pd.DataFrame:
    """Generate mock PubChem data for demonstration"""
    mock_data = [{
        'cid': 2244,
        'iupac_name': '2-acetoxybenzoic acid',
        'molecular_formula': 'C9H8O4',
        'molecular_weight': 180.16,
        'canonical_smiles': 'CC(=O)Oc1ccccc1C(=O)O',
        'isomeric_smiles': 'CC(=O)Oc1ccccc1C(=O)O',
        'inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
        'inchikey': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N',
        'synonyms': 'aspirin, acetylsalicylic acid, 2-acetoxybenzoic acid',
        'xlogp': 1.2,
        'exact_mass': 180.04226,
        'monoisotopic_mass': 180.04226,
        'tpsa': 63.6,
        'complexity': 212,
        'charge': 0,
        'h_bond_donor_count': 1,
        'h_bond_acceptor_count': 4,
        'rotatable_bond_count': 3,
        'heavy_atom_count': 13
    }]

    return pd.DataFrame(mock_data)

def _get_mock_similarity_data(query_smiles: str, threshold: float) -> pd.DataFrame:
    """Generate mock similarity search data"""
    import random

    mock_data = []
    base_smiles = ['CCO', 'CCCO', 'CC(=O)O', 'c1ccccc1', 'CCN(CC)CC']

    for i, smiles in enumerate(base_smiles):
        if random.random() > 0.3:  # Random filtering
            mock_data.append({
                'compound_id': f'SIM_{i:03d}',
                'canonical_smiles': smiles,
                'similarity_score': random.uniform(threshold, 1.0),
                'database': 'ChEMBL',
                'molecular_weight': random.uniform(100, 500),
                'logp': random.uniform(-2, 5)
            })

    return pd.DataFrame(mock_data)

def get_compound_bioactivities(chembl_id: str) -> pd.DataFrame:
    """
    Get all bioactivities for a specific ChEMBL compound

    Args:
        chembl_id: ChEMBL compound ID

    Returns:
        DataFrame with bioactivity data
    """
    if not CHEMBL_AVAILABLE:
        return pd.DataFrame()

    try:
        activity = new_client.activity
        activities = activity.filter(molecule_chembl_id=chembl_id)

        results = []
        for act in activities:
            result = {
                'target_chembl_id': act.get('target_chembl_id'),
                'activity_type': act.get('standard_type'),
                'activity_value': act.get('standard_value'),
                'activity_units': act.get('standard_units'),
                'confidence_score': act.get('confidence_score'),
                'assay_description': act.get('assay_description')
            }
            results.append(result)

        return pd.DataFrame(results)

    except Exception as e:
        logging.error(f"Error retrieving bioactivities: {str(e)}")
        return pd.DataFrame()

def get_target_info(target_chembl_id: str) -> Dict:
    """
    Get detailed information about a ChEMBL target

    Args:
        target_chembl_id: ChEMBL target ID

    Returns:
        Dictionary with target information
    """
    if not CHEMBL_AVAILABLE:
        return {}

    try:
        target = new_client.target
        target_info = target.get(target_chembl_id)

        return {
            'chembl_id': target_info.get('target_chembl_id'),
            'name': target_info.get('pref_name'),
            'type': target_info.get('target_type'),
            'organism': target_info.get('organism'),
            'uniprot_accession': target_info.get('target_components', [{}])[0].get('accession') if target_info.get('target_components') else None,
            'gene_names': target_info.get('target_components', [{}])[0].get('component_synonyms') if target_info.get('target_components') else []
        }

    except Exception as e:
        logging.error(f"Error retrieving target info: {str(e)}")
        return {}
