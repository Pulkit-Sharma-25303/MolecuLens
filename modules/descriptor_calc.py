"""
Molecular descriptor calculation utilities
"""
from rdkit import Chem
from mordred import Calculator, descriptors
import pandas as pd

def get_mordred_df(smiles_list):
    """
    Takes a list of SMILES strings and returns a Pandas DataFrame of Mordred descriptors.
    """
    calc = Calculator(descriptors, ignore_3D=True)
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    desc = calc.pandas(mols)
    desc = desc.select_dtypes(include=['number'])  # Only numeric columns
    desc.reset_index(drop=True, inplace=True)
    return desc

import pandas as pd
import numpy as np
import logging
from typing import List

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

try:
    from mordred import Calculator, descriptors
    MORDRED_AVAILABLE = True
except ImportError:
    MORDRED_AVAILABLE = False

from .molecular_utils import smiles_to_mol

def calculate_basic_properties(smiles_list: List[str]) -> pd.DataFrame:
    """Calculate basic molecular properties for a list of SMILES"""
    if not RDKIT_AVAILABLE:
        logging.warning("RDKit not available")
        return pd.DataFrame()

    results = []
    for smiles in smiles_list:
        mol = smiles_to_mol(smiles)
        if mol is None:
            continue
        try:
            properties = {
                'SMILES': smiles,
                'MW': Descriptors.MolWt(mol),
                'LogP': Crippen.MolLogP(mol),
                'HBD': rdMolDescriptors.CalcNumHBD(mol),
                'HBA': rdMolDescriptors.CalcNumHBA(mol),
                'TPSA': rdMolDescriptors.CalcTPSA(mol),
                'RotBonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
                'AromaticRings': rdMolDescriptors.CalcNumAromaticRings(mol),
                'SaturatedRings': rdMolDescriptors.CalcNumSaturatedRings(mol),
                'HeavyAtoms': mol.GetNumHeavyAtoms(),
                'FormalCharge': Chem.GetFormalCharge(mol)
            }
            results.append(properties)
        except Exception as e:
            logging.error(f"Error calculating properties for {smiles}: {str(e)}")
    return pd.DataFrame(results)

def calculate_lipinski(smiles_list: List[str]) -> pd.DataFrame:
    """Calculate Lipinski Rule of Five parameters."""
    if not RDKIT_AVAILABLE:
        return pd.DataFrame()
    results = []
    for smiles in smiles_list:
        mol = smiles_to_mol(smiles)
        if mol is None:
            continue
        try:
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            violations = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)
            result = {
                'SMILES': smiles,
                'MW': mw,
                'LogP': logp,
                'HBD': hbd,
                'HBA': hba,
                'MW_Pass': mw <= 500,
                'LogP_Pass': logp <= 5,
                'HBD_Pass': hbd <= 5,
                'HBA_Pass': hba <= 10,
                'Lipinski_Violations': violations,
                'Drug_Like': violations <= 1
            }
            results.append(result)
        except Exception as e:
            logging.error(f"Error calculating Lipinski for {smiles}: {str(e)}")
    return pd.DataFrame(results)

def calculate_admet_properties(smiles_list: List[str]) -> pd.DataFrame:
    """Calculate ADMET-related properties (with mock/simple predictions for demo)."""
    if not RDKIT_AVAILABLE:
        return pd.DataFrame()
    results = []
    for smiles in smiles_list:
        mol = smiles_to_mol(smiles)
        if mol is None:
            continue
        try:
            MW = Descriptors.MolWt(mol)
            TPSA = rdMolDescriptors.CalcTPSA(mol)
            LogP = Crippen.MolLogP(mol)
            properties = {
                'SMILES': smiles,
                'MW': MW,
                'LogP': LogP,
                'TPSA': TPSA,
                'HBD': rdMolDescriptors.CalcNumHBD(mol),
                'HBA': rdMolDescriptors.CalcNumHBA(mol),
                'RotBonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
                'Absorption_Pred': 'High' if TPSA < 140 and MW < 500 else 'Low',
                'BBB_Pred': 'Yes' if TPSA < 90 and MW < 450 else 'No',
                'Solubility_Pred': 'Soluble' if LogP < 3 else 'Poorly Soluble',
                'CYP2D6_Substrate': np.random.choice(['Yes', 'No']),
                'hERG_Risk': 'High' if LogP > 4 and MW > 300 else 'Low'
            }
            results.append(properties)
        except Exception as e:
            logging.error(f"Error calculating ADMET for {smiles}: {str(e)}")
    return pd.DataFrame(results)

def calculate_descriptors(smiles_list: List[str]) -> pd.DataFrame:
    """Calculate comprehensive descriptors using Mordred, fallback to RDKit if unavailable."""
    if not MORDRED_AVAILABLE:
        logging.warning("Mordred not available, using RDKit descriptors only")
        return calculate_rdkit_descriptors(smiles_list)
    try:
        mols = [smiles_to_mol(smiles) for smiles in smiles_list]
        valid = [(mol, smi) for mol, smi in zip(mols, smiles_list) if mol is not None]
        if not valid:
            return pd.DataFrame()
        valid_mols, valid_smiles = zip(*valid)
        calc = Calculator(descriptors, ignore_3D=True)
        desc_df = calc.pandas(valid_mols)
        desc_df.insert(0, 'SMILES', valid_smiles)
        desc_df = desc_df.dropna(axis=1, how='all')
        return desc_df
    except Exception as e:
        logging.error(f"Error calculating Mordred descriptors: {str(e)}")
        return calculate_rdkit_descriptors(smiles_list)

def calculate_rdkit_descriptors(smiles_list: List[str]) -> pd.DataFrame:
    """Fallback: Calculate all available RDKit descriptors."""
    if not RDKIT_AVAILABLE:
        return pd.DataFrame()
    results = []
    descriptor_names = [desc[0] for desc in Descriptors._descList]
    for smiles in smiles_list:
        mol = smiles_to_mol(smiles)
        if mol is None:
            continue
        try:
            props = {'SMILES': smiles}
            for desc_name in descriptor_names:
                try:
                    func = getattr(Descriptors, desc_name)
                    props[desc_name] = func(mol)
                except Exception:
                    props[desc_name] = np.nan
            results.append(props)
        except Exception as e:
            logging.error(f"Error calculating RDKit descriptors for {smiles}: {str(e)}")
    return pd.DataFrame(results)

def get_descriptor_names() -> List[str]:
    """List all available descriptor names (RDKit and Mordred if installed)."""
    names = []
    if RDKIT_AVAILABLE:
        names.extend([desc[0] for desc in Descriptors._descList])
    if MORDRED_AVAILABLE:
        calc = Calculator(descriptors, ignore_3D=True)
        names.extend([str(desc) for desc in calc.descriptors])
    return sorted(set(names))

def filter_descriptors(df: pd.DataFrame, remove_constant: bool = True,
                       remove_correlated: bool = True,
                       correlation_threshold: float = 0.95) -> pd.DataFrame:
    """Filter constant and highly correlated descriptor columns."""
    filtered_df = df.copy()
    numeric_cols = filtered_df.select_dtypes(include=[np.number]).columns
    if remove_constant:
        constant_cols = [col for col in numeric_cols if filtered_df[col].nunique() <= 1]
        filtered_df = filtered_df.drop(columns=constant_cols)
        logging.info(f"Removed {len(constant_cols)} constant descriptors")
    if remove_correlated:
        numeric_cols = filtered_df.select_dtypes(include=[np.number]).columns
        corr_matrix = filtered_df[numeric_cols].corr().abs()
        upper_tri = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
        high_corr_cols = [
            column for column in upper_tri.columns
            if any(upper_tri[column] > correlation_threshold)
        ]
        filtered_df = filtered_df.drop(columns=high_corr_cols)
        logging.info(f"Removed {len(high_corr_cols)} highly correlated descriptors")
    return filtered_df

def calculate_molecular_fingerprints(smiles_list: List[str], fp_type: str = 'morgan',
                                     radius: int = 2, n_bits: int = 2048) -> pd.DataFrame:
    """Calculate molecular fingerprints for input SMILES."""
    if not RDKIT_AVAILABLE:
        return pd.DataFrame()
    try:
        from rdkit.Chem import rdFingerprintGenerator
        results = []
        for smiles in smiles_list:
            mol = smiles_to_mol(smiles)
            if mol is None:
                continue
            if fp_type == 'morgan':
                fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
            elif fp_type == 'rdkit':
                fp_gen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=n_bits)
            else:
                fp_gen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=n_bits)
            fp = fp_gen.GetFingerprint(mol)
            fp_array = np.array(fp)
            fp_dict = {'SMILES': smiles}
            fp_dict.update({f'Bit_{i}': int(bit) for i, bit in enumerate(fp_array)})
            results.append(fp_dict)
        return pd.DataFrame(results)
    except Exception as e:
        logging.error(f"Error in fingerprint calculation: {str(e)}")
        return pd.DataFrame()
