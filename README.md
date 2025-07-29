# MolecuLens - AI-Powered Drug Discovery Platform

🧬 **A comprehensive web application for computational drug discovery built entirely with Python and Streamlit**

## Overview

MolecuLens is a full-featured drug discovery platform that integrates molecular visualization, property calculation, database search, virtual screening, and molecular dynamics simulation into a single, user-friendly web interface. Built entirely with Python and the scientific Python ecosystem, it provides researchers with powerful tools for computational chemistry and drug design.

## Features

### 🧪 Molecular Editor
- Interactive molecular structure drawing and editing
- SMILES input/output with validation
- 3D molecular visualization using Stmol and py3DMol
- Real-time property calculation
- Support for multiple input formats (SMILES, MOL, SDF)

### ⚗️ Properties Calculator
- Comprehensive molecular descriptor calculation (1600+ descriptors via Mordred)
- Lipinski's Rule of Five analysis
- ADMET property predictions
- Batch processing capabilities
- Interactive property visualization
- Export results in multiple formats

### 🔍 Database Search
- ChEMBL bioactivity database integration
- PubChem chemical information queries
- Target-based and compound-based searches
- Similarity search with configurable thresholds
- Substructure filtering
- Export search results

### 📊 Virtual Screening
- Compound library screening with customizable filters
- Lipinski Rule of Five filtering
- PAINS (Pan-Assay Interference Compounds) detection
- Substructure and similarity-based filtering
- Property-based filtering (MW, LogP, etc.)
- Interactive results visualization and ranking

### 🎯 Molecular Dynamics Simulation
- OpenMM integration for MD simulations
- Multiple force fields (AMBER, CHARMM, OPLS-AA)
- Trajectory analysis with MDAnalysis
- RMSD, RMSF, and other property calculations
- GPU acceleration support
- Interactive result visualization

## Technology Stack

### Core Framework
- **Streamlit** - Web application framework
- **Python 3.8+** - Primary programming language

### Cheminformatics
- **RDKit** - Molecular handling and property calculation
- **Mordred** - Comprehensive descriptor calculation
- **OpenEye** - Advanced cheminformatics tools (optional)

### Visualization
- **Stmol** - 3D molecular visualization in Streamlit
- **py3DMol** - WebGL-based molecular rendering
- **mols2grid** - Interactive molecular grids
- **Plotly** - Interactive plots and charts

### Molecular Simulation
- **OpenMM** - Molecular dynamics engine
- **MDAnalysis** - Trajectory analysis
- **NumPy/SciPy** - Numerical computations

### Database Integration
- **ChEMBL Web Resource Client** - ChEMBL database access
- **PubChemPy** - PubChem database queries
- **Pandas** - Data manipulation and analysis

## Installation

### Prerequisites
- Python 3.8 or later
- pip package manager

### Quick Start

1. **Clone or download the project:**
   ```bash
   # If you have the zip file, extract it
   unzip moleculens.zip
   cd moleculens
   ```

2. **Create a virtual environment (recommended):**
   ```bash
   python -m venv venv

   # On Linux/Mac:
   source venv/bin/activate

   # On Windows:
   venv\Scripts\activate
   ```

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

4. **Run the application:**
   ```bash
   streamlit run app.py
   ```

5. **Open in browser:**
   The application will automatically open in your default browser at `http://localhost:8501`

### Alternative Installation Methods

#### Using Conda
```bash
conda create -n moleculens python=3.9
conda activate moleculens
pip install -r requirements.txt
streamlit run app.py
```

#### Docker (Advanced)
```bash
# Build Docker image
docker build -t moleculens .

# Run container
docker run -p 8501:8501 moleculens
```

## Usage

### Getting Started
1. **Home Page**: Navigate through the application using the sidebar menu
2. **Molecular Editor**: Start by drawing or inputting a molecule structure
3. **Properties**: Calculate molecular descriptors and drug-like properties
4. **Database Search**: Search ChEMBL and PubChem for bioactivity data
5. **Virtual Screening**: Screen compound libraries with customizable filters
6. **MD Simulation**: Set up and run molecular dynamics simulations

### Example Workflows

#### Drug-Like Property Analysis
1. Go to **Molecular Editor**
2. Input SMILES: `CCc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc4nccc(-c5cccnc5)n4`
3. Navigate to **Properties Calculator**
4. Select "Lipinski Rule of Five" analysis
5. Review results and export data

#### Target-Based Database Search
1. Go to **Database Search**
2. Select "Target Search"
3. Enter target: "EGFR"
4. Set activity filters (e.g., IC50 < 100 nM)
5. Export results for further analysis

#### Virtual Library Screening
1. Go to **Virtual Screening**
2. Upload compound library or use sample data
3. Configure filters (Lipinski, PAINS, similarity)
4. Run screening and analyze results
5. Export hit compounds

## Configuration

### Environment Variables
Create a `.env` file in the project root for configuration:
```
CHEMBL_BASE_URL=https://www.ebi.ac.uk/chembl/api/data
PUBCHEM_BASE_URL=https://pubchem.ncbi.nlm.nih.gov/rest/pug
OPENMM_CUDA_PRECISION=mixed
```

### Custom Settings
Modify `config.py` for advanced configuration:
- Database connection settings
- Simulation parameters
- Visualization preferences
- Export formats

## Development

### Project Structure
```
moleculens/
├── app.py                      # Main Streamlit application
├── requirements.txt            # Python dependencies
├── pages/                      # Streamlit pages
│   ├── 1_🏠_Home.py           # Home page
│   ├── 2_🧪_Molecular_Editor.py # Molecular editor
│   ├── 3_⚗️_Properties.py      # Properties calculator
│   ├── 4_🔍_Database_Search.py # Database search
│   ├── 5_📊_Virtual_Screening.py # Virtual screening
│   └── 6_🎯_MD_Simulation.py   # MD simulation
├── modules/                    # Core functionality modules
│   ├── __init__.py
│   ├── molecular_utils.py      # Molecular handling utilities
│   ├── descriptor_calc.py      # Descriptor calculations
│   ├── database_client.py      # Database clients
│   ├── visualization.py        # Visualization utilities
│   ├── screening_utils.py      # Screening algorithms
│   └── md_simulator.py         # MD simulation wrapper
└── assets/                     # Static assets
    └── logo.png               # Application logo
```

### Adding New Features
1. Create new modules in the `modules/` directory
2. Add new pages in the `pages/` directory
3. Update `requirements.txt` for new dependencies
4. Follow existing code patterns and documentation

### Testing
```bash
# Install testing dependencies
pip install pytest pytest-cov

# Run tests
pytest tests/

# Generate coverage report
pytest --cov=modules tests/
```

## Deployment

### Streamlit Community Cloud (Recommended)
1. Push code to GitHub repository
2. Connect to Streamlit Community Cloud
3. Deploy with one click
4. Share the public URL

### Other Deployment Options
- **Heroku**: Use `Procfile` for Heroku deployment
- **AWS/GCP/Azure**: Deploy using container services
- **Local Server**: Run on local network for team access

## API Reference

### Core Modules

#### molecular_utils
```python
from modules.molecular_utils import smiles_to_mol, validate_smiles

mol = smiles_to_mol("CCO")
is_valid, message = validate_smiles("CCO")
```

#### descriptor_calc
```python
from modules.descriptor_calc import calculate_lipinski

results = calculate_lipinski(["CCO", "CCCO"])
```

#### database_client
```python
from modules.database_client import search_chembl_by_target

results = search_chembl_by_target("EGFR", limit=100)
```

## Contributing

We welcome contributions!

### Development Setup
1. Fork the repository
2. Create a feature branch
3. Make changes and add tests
4. Submit a pull request

### Reporting Issues
- Use the GitHub issue tracker
- Provide detailed error messages
- Include system information
- Describe steps to reproduce


## Acknowledgments

- **RDKit** - Open-source cheminformatics toolkit
- **OpenMM** - Molecular simulation engine
- **Streamlit** - Web application framework
- **ChEMBL** - Medicinal chemistry database
- **PubChem** - Chemical information database

## Changelog

### Version 1.0.0 (2025)
- Initial release
- Complete drug discovery pipeline
- Six main modules with comprehensive functionality
- Documentation and examples

---

**MolecuLens** - Advancing drug discovery through computational chemistry and AI
