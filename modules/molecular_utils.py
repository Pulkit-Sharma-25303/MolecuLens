"""
Molecular utilities for handling chemical structures
"""

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

import logging

def smiles_to_mol(smiles):
    """
    Convert SMILES string to RDKit molecule object

    Args:
        smiles (str): SMILES string

    Returns:
        rdkit.Chem.Mol: RDKit molecule object or None if invalid
    """
    if not RDKIT_AVAILABLE:
        logging.warning("RDKit not available. Install rdkit-pypi for full functionality.")
        return None

    if not smiles or not isinstance(smiles, str):
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Sanitize molecule
        Chem.SanitizeMol(mol)
        return mol
    except Exception as e:
        logging.error(f"Error converting SMILES {smiles}: {str(e)}")
        return None

def mol_to_smiles(mol):
    """
    Convert RDKit molecule to canonical SMILES

    Args:
        mol: RDKit molecule object

    Returns:
        str: Canonical SMILES string
    """
    if not RDKIT_AVAILABLE or mol is None:
        return ""

    try:
        return Chem.MolToSmiles(mol)
    except Exception as e:
        logging.error(f"Error converting molecule to SMILES: {str(e)}")
        return ""

def validate_smiles(smiles):
    """
    Validate SMILES string

    Args:
        smiles (str): SMILES string to validate

    Returns:
        tuple: (is_valid, error_message)
    """
    if not smiles or not isinstance(smiles, str):
        return False, "Empty or invalid SMILES string"

    if not RDKIT_AVAILABLE:
        return True, "RDKit not available for validation"

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES syntax"

        # Check for common issues
        if mol.GetNumAtoms() == 0:
            return False, "No atoms in molecule"

        # Try to sanitize
        Chem.SanitizeMol(mol)
        return True, "Valid SMILES"

    except Exception as e:
        return False, f"Validation error: {str(e)}"

def validate_smiles_list(smiles_list):
    """
    Validate list of SMILES strings

    Args:
        smiles_list (list): List of SMILES strings

    Returns:
        tuple: (valid_smiles, invalid_smiles)
    """
    valid_smiles = []
    invalid_smiles = []

    for smiles in smiles_list:
        is_valid, _ = validate_smiles(smiles)
        if is_valid:
            valid_smiles.append(smiles)
        else:
            invalid_smiles.append(smiles)

    return valid_smiles, invalid_smiles

def calculate_molecular_weight(mol):
    """
    Calculate molecular weight

    Args:
        mol: RDKit molecule object

    Returns:
        float: Molecular weight in Da
    """
    if not RDKIT_AVAILABLE or mol is None:
        return 0.0

    try:
        return Descriptors.MolWt(mol)
    except Exception:
        return 0.0

def calculate_logp(mol):
    """
    Calculate LogP (partition coefficient)

    Args:
        mol: RDKit molecule object

    Returns:
        float: LogP value
    """
    if not RDKIT_AVAILABLE or mol is None:
        return 0.0

    try:
        return Crippen.MolLogP(mol)
    except Exception:
        return 0.0

def calculate_hbd_hba(mol):
    """
    Calculate hydrogen bond donors and acceptors

    Args:
        mol: RDKit molecule object

    Returns:
        tuple: (donors, acceptors)
    """
    if not RDKIT_AVAILABLE or mol is None:
        return 0, 0

    try:
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        return hbd, hba
    except Exception:
        return 0, 0

def calculate_tpsa(mol):
    """
    Calculate topological polar surface area

    Args:
        mol: RDKit molecule object

    Returns:
        float: TPSA in Ų
    """
    if not RDKIT_AVAILABLE or mol is None:
        return 0.0

    try:
        return rdMolDescriptors.CalcTPSA(mol)
    except Exception:
        return 0.0

def calculate_rotatable_bonds(mol):
    """
    Calculate number of rotatable bonds

    Args:
        mol: RDKit molecule object

    Returns:
        int: Number of rotatable bonds
    """
    if not RDKIT_AVAILABLE or mol is None:
        return 0

    try:
        return rdMolDescriptors.CalcNumRotatableBonds(mol)
    except Exception:
        return 0

def add_hydrogens(mol):
    """
    Add explicit hydrogens to molecule

    Args:
        mol: RDKit molecule object

    Returns:
        rdkit.Chem.Mol: Molecule with explicit hydrogens
    """
    if not RDKIT_AVAILABLE or mol is None:
        return mol

    try:
        return Chem.AddHs(mol)
    except Exception:
        return mol

def remove_hydrogens(mol):
    """
    Remove explicit hydrogens from molecule

    Args:
        mol: RDKit molecule object

    Returns:
        rdkit.Chem.Mol: Molecule without explicit hydrogens
    """
    if not RDKIT_AVAILABLE or mol is None:
        return mol

    try:
        return Chem.RemoveHs(mol)
    except Exception:
        return mol

def generate_conformer(mol, num_confs=1):
    """
    Generate 3D conformer for molecule

    Args:
        mol: RDKit molecule object
        num_confs (int): Number of conformers to generate

    Returns:
        rdkit.Chem.Mol: Molecule with 3D coordinates
    """
    if not RDKIT_AVAILABLE or mol is None:
        return mol

    try:
        from rdkit.Chem import AllChem

        # Add hydrogens and generate conformer
        mol_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_h, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol_h)

        return mol_h
    except Exception as e:
        logging.error(f"Error generating conformer: {str(e)}")
        return mol

def calculate_basic_properties(mol):
    """
    Calculate basic molecular properties

    Args:
        mol: RDKit molecule object

    Returns:
        dict: Dictionary of basic properties
    """
    if not RDKIT_AVAILABLE or mol is None:
        return {}

    try:
        mw = calculate_molecular_weight(mol)
        logp = calculate_logp(mol)
        hbd, hba = calculate_hbd_hba(mol)
        tpsa = calculate_tpsa(mol)
        rotbonds = calculate_rotatable_bonds(mol)

        return {
            'MW': mw,
            'LogP': logp,
            'HBD': hbd,
            'HBA': hba,
            'TPSA': tpsa,
            'RotBonds': rotbonds
        }
    except Exception as e:
        logging.error(f"Error calculating properties: {str(e)}")
        return {}

# Utility functions for file handling
def read_sdf_file(file_path):
    """
    Read molecules from SDF file

    Args:
        file_path (str): Path to SDF file

    Returns:
        list: List of RDKit molecule objects
    """
    if not RDKIT_AVAILABLE:
        return []

    try:
        supplier = Chem.SDMolSupplier(file_path)
        molecules = [mol for mol in supplier if mol is not None]
        return molecules
    except Exception as e:
        logging.error(f"Error reading SDF file: {str(e)}")
        return []

def write_sdf_file(molecules, file_path):
    """
    Write molecules to SDF file

    Args:
        molecules (list): List of RDKit molecule objects
        file_path (str): Output file path

    Returns:
        bool: Success status
    """
    if not RDKIT_AVAILABLE:
        return False

    try:
        writer = Chem.SDWriter(file_path)
        for mol in molecules:
            if mol is not None:
                writer.write(mol)
        writer.close()
        return True
    except Exception as e:
        logging.error(f"Error writing SDF file: {str(e)}")
        return False
