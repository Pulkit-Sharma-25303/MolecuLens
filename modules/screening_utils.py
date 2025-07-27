"""
Screening utilities for virtual screening and molecular filtering
"""

import pandas as pd
import numpy as np
import logging
from typing import List, Optional

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, DataStructs
    from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

from .molecular_utils import smiles_to_mol

def lipinski_flags(smiles: str) -> bool:
    """
    Check if a molecule passes Lipinski's Rule of Five.
    
    Args:
        smiles: SMILES string
        
    Returns:
        bool: True if molecule passes Lipinski rules (≤1 violation), False otherwise
    """
    if not RDKIT_AVAILABLE:
        return False
        
    mol = smiles_to_mol(smiles)
    if mol is None:
        return False
        
    try:
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(NumHBA(mol))
        
        violations = 0
        violations += 1 if mw > 500 else 0
        violations += 1 if logp > 5 else 0
        violations += 1 if hbd > 5 else 0
        violations += 1 if hba > 10 else 0
        
        return violations <= 1
        
    except Exception as e:
        logging.error(f"Error calculating Lipinski for {smiles}: {str(e)}")
        return False

def filter_by_lipinski(df: pd.DataFrame, smiles_column: str = 'SMILES') -> pd.DataFrame:
    """
    Filter DataFrame by Lipinski Rule of Five.
    
    Args:
        df: DataFrame containing SMILES
        smiles_column: Name of column containing SMILES strings
        
    Returns:
        DataFrame with Lipinski-compliant molecules
    """
    if smiles_column not in df.columns:
        return df
        
    df['Lipinski_Pass'] = df[smiles_column].apply(lipinski_flags)
    return df[df['Lipinski_Pass']].copy()

def calculate_lipinski_properties(smiles: str) -> dict:
    """
    Calculate detailed Lipinski properties for a single molecule.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dictionary with Lipinski properties
    """
    if not RDKIT_AVAILABLE:
        return {}
        
    mol = smiles_to_mol(smiles)
    if mol is None:
        return {}
        
    try:
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        
        violations = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)
        
        return {
            'MW': mw,
            'LogP': logp,
            'HBD': hbd,
            'HBA': hba,
            'MW_Pass': mw <= 500,
            'LogP_Pass': logp <= 5,
            'HBD_Pass': hbd <= 5,
            'HBA_Pass': hba <= 10,
            'Violations': violations,
            'Drug_Like': violations <= 1
        }
        
    except Exception as e:
        logging.error(f"Error calculating Lipinski properties for {smiles}: {str(e)}")
        return {}

def substructure_filter(df: pd.DataFrame, smarts_pattern: str, 
                       smiles_column: str = 'SMILES') -> pd.DataFrame:
    """
    Filter DataFrame by substructure match.
    
    Args:
        df: DataFrame containing SMILES
        smarts_pattern: SMARTS pattern to match
        smiles_column: Name of column containing SMILES strings
        
    Returns:
        DataFrame with molecules containing the substructure
    """
    if not RDKIT_AVAILABLE or smiles_column not in df.columns:
        return df
        
    try:
        pattern_mol = Chem.MolFromSmarts(smarts_pattern)
        if pattern_mol is None:
            logging.warning(f"Invalid SMARTS pattern: {smarts_pattern}")
            return df
            
        def has_substructure(smiles):
            mol = smiles_to_mol(smiles)
            if mol is None:
                return False
            return mol.HasSubstructMatch(pattern_mol)
            
        df['Has_Substructure'] = df[smiles_column].apply(has_substructure)
        return df[df['Has_Substructure']].copy()
        
    except Exception as e:
        logging.error(f"Error in substructure filtering: {str(e)}")
        return df

def similarity_screen(df: pd.DataFrame, reference_smiles: str, 
                     threshold: float = 0.7, smiles_column: str = 'SMILES') -> pd.DataFrame:
    """
    Filter DataFrame by Tanimoto similarity to reference molecule.
    
    Args:
        df: DataFrame containing SMILES
        reference_smiles: Reference SMILES for similarity comparison
        threshold: Minimum similarity threshold (0.0 to 1.0)
        smiles_column: Name of column containing SMILES strings
        
    Returns:
        DataFrame with molecules above similarity threshold
    """
    if not RDKIT_AVAILABLE or smiles_column not in df.columns:
        return df
        
    try:
        ref_mol = smiles_to_mol(reference_smiles)
        if ref_mol is None:
            logging.warning(f"Invalid reference SMILES: {reference_smiles}")
            return df
            
        ref_fp = GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=2048)
        
        def calculate_similarity(smiles):
            mol = smiles_to_mol(smiles)
            if mol is None:
                return 0.0
            fp = GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            return DataStructs.TanimotoSimilarity(ref_fp, fp)
            
        df['Similarity'] = df[smiles_column].apply(calculate_similarity)
        return df[df['Similarity'] >= threshold].copy()
        
    except Exception as e:
        logging.error(f"Error in similarity screening: {str(e)}")
        return df

def filter_by_properties(df: pd.DataFrame, smiles_column: str = 'SMILES',
                        mw_range: tuple = (0, 1000), 
                        logp_range: tuple = (-10, 10),
                        hbd_range: tuple = (0, 20),
                        hba_range: tuple = (0, 30)) -> pd.DataFrame:
    """
    Filter DataFrame by multiple molecular properties.
    
    Args:
        df: DataFrame containing SMILES
        smiles_column: Name of column containing SMILES strings
        mw_range: Molecular weight range (min, max)
        logp_range: LogP range (min, max)
        hbd_range: H-bond donor range (min, max)
        hba_range: H-bond acceptor range (min, max)
        
    Returns:
        Filtered DataFrame
    """
    if not RDKIT_AVAILABLE or smiles_column not in df.columns:
        return df
        
    def passes_filters(smiles):
        mol = smiles_to_mol(smiles)
        if mol is None:
            return False
            
        try:
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            
            return (mw_range[0] <= mw <= mw_range[1] and
                   logp_range[0] <= logp <= logp_range[1] and
                   hbd_range[0] <= hbd <= hbd_range[1] and
                   hba_range[0] <= hba <= hba_range[1])
                   
        except Exception:
            return False
            
    df['Passes_Filters'] = df[smiles_column].apply(passes_filters)
    return df[df['Passes_Filters']].copy()

def lead_like_filter(smiles: str) -> bool:
    """
    Check if molecule passes lead-like criteria.
    Lead-like: MW 250-350, LogP 1-3, HBD ≤3, HBA ≤6, RotBonds ≤7
    
    Args:
        smiles: SMILES string
        
    Returns:
        bool: True if lead-like, False otherwise
    """
    if not RDKIT_AVAILABLE:
        return False
        
    mol = smiles_to_mol(smiles)
    if mol is None:
        return False
        
    try:
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rotbonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        
        return (250 <= mw <= 350 and
               1 <= logp <= 3 and
               hbd <= 3 and
               hba <= 6 and
               rotbonds <= 7)
               
    except Exception as e:
        logging.error(f"Error calculating lead-like properties for {smiles}: {str(e)}")
        return False

def fragment_like_filter(smiles: str) -> bool:
    """
    Check if molecule passes fragment-like criteria.
    Fragment-like: MW ≤300, LogP ≤3, HBD ≤3, HBA ≤6, RotBonds ≤3
    
    Args:
        smiles: SMILES string
        
    Returns:
        bool: True if fragment-like, False otherwise
    """
    if not RDKIT_AVAILABLE:
        return False
        
    mol = smiles_to_mol(smiles)
    if mol is None:
        return False
        
    try:
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rotbonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        
        return (mw <= 300 and
               logp <= 3 and
               hbd <= 3 and
               hba <= 6 and
               rotbonds <= 3)
               
    except Exception as e:
        logging.error(f"Error calculating fragment-like properties for {smiles}: {str(e)}")
        return False

def pains_filter(smiles: str) -> bool:
    """
    Simple PAINS (Pan-Assay Interference Compounds) filter.
    This is a simplified version - in practice, you'd want a more comprehensive filter.
    
    Args:
        smiles: SMILES string
        
    Returns:
        bool: True if passes PAINS filter, False if flagged as PAINS
    """
    if not RDKIT_AVAILABLE:
        return True
        
    mol = smiles_to_mol(smiles)
    if mol is None:
        return False
        
    # Simple PAINS patterns (this is a very basic implementation)
    pains_patterns = [
        'c1ccc2c(c1)c(=O)c(=O)[nH]2',  # Quinone
        'N=Nc1ccccc1',  # Azo compound
        'c1ccc(cc1)S(=O)(=O)N',  # Sulfonamide (basic pattern)
    ]
    
    try:
        for pattern in pains_patterns:
            pattern_mol = Chem.MolFromSmarts(pattern)
            if pattern_mol and mol.HasSubstructMatch(pattern_mol):
                return False
        return True
        
    except Exception as e:
        logging.error(f"Error in PAINS filtering for {smiles}: {str(e)}")
        return True

def calculate_diversity_score(smiles_list: List[str]) -> float:
    """
    Calculate diversity score for a list of molecules using Morgan fingerprints.
    
    Args:
        smiles_list: List of SMILES strings
        
    Returns:
        float: Average pairwise Tanimoto distance (0 = identical, 1 = completely diverse)
    """
    if not RDKIT_AVAILABLE or len(smiles_list) < 2:
        return 0.0
        
    try:
        fps = []
        for smiles in smiles_list:
            mol = smiles_to_mol(smiles)
            if mol is not None:
                fp = GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                fps.append(fp)
                
        if len(fps) < 2:
            return 0.0
            
        distances = []
        for i in range(len(fps)):
            for j in range(i+1, len(fps)):
                similarity = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                distance = 1.0 - similarity
                distances.append(distance)
                
        return np.mean(distances) if distances else 0.0
        
    except Exception as e:
        logging.error(f"Error calculating diversity score: {str(e)}")
        return 0.0
