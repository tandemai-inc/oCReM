# oCReM

This repository provides the database and tools for the paper "Expanding accessible chemical space for fragment-based enumeration by orders of magnitude through optimization of the CReM framework". 

It contains the data of modified CrEM library and the original ChEMBL fragment library, along with tools for molecular fragmentation, fragment database management, and structure generation.

## Authorlist

- Shaojin Hu, Qinyu Chen, Yinhui Yi, Paul Pilot, James Xu, Abir Ganguly and Albert C. Pan

## Data Access
Original ChemBL fragment library: [http://www.qsar4u.com/pages/crem.php](https://www.qsar4u.com/pages/crem.php)

## Installation

### Create conda environment
Make sure you have installed miniconda.
```bash
conda env create -f environment.yml
conda activate crem
```

## Overview

oCReM is a framework for fragment-based molecular design that includes:

- **Molecular Fragmentation**: Breaking down molecules into structural fragments
- **Fragment Database Management**: Storing and retrieving fragments in SQLite or PostgreSQL databases
- **Structure Generation**: Generating new molecules through fragment-based assembly
  - Mutate: Replace fragments in a molecule
  - Grow: Extend a molecule with new fragments
  - Link: Connect multiple molecules using linker fragments

## Directory Structure

```
oCReM/
├── example/              # Example Jupyter notebooks
│   ├── fragmentation_to_db.ipynb  # Fragmentation and database creation example
│   └── structure_generation.ipynb # Structure generation example
├── ta_gen/               # Core functionality
│   ├── bin/              # Command-line tools
│   ├── crem/             # Core CReM implementation
│   ├── db/               # Database management
│   └── utils/            # Utility functions
├── LICENSE               # License file
├── README.md             # This README
└── environment.yml       # Conda environment configuration
```

## Usage

### 1. Molecular Fragmentation

Fragment molecules and save results to a file:

```bash
python ta_gen/bin/fragmentation.py --input test.smi --out test_frag.csv --mode 0 --ncpu 10 --radius 3
```

Fragment molecules and save results to a SQLite database:

```bash
python ta_gen/bin/fragmentation.py --input test.smi --mode 0 --ncpu 10 --radius 3 --use_db --db_type sqlite --db_path test.db
```

Fragment molecules and save results to a PostgreSQL database:

```bash
# First create a PostgreSQL configuration file (e.g., test.ini)
python ta_gen/bin/fragmentation.py --input test.smi --mode 0 --ncpu 10 --radius 3 --use_db --db_type postgres --ini_file test.ini
```

For detailed examples and step-by-step instructions, please refer to the Jupyter notebooks in the `example/` directory, especially `fragmentation_to_db.ipynb` which covers the complete fragmentation process and database storage options.

### 2. Structure Generation

#### Mutate Molecule

```python
from rdkit import Chem
from ta_gen.crem.crem import mutate_mol
from ta_gen.db import create_db_manager

m = Chem.MolFromSmiles('c1cc(OC)ccc1C')  # toluene
db_manager = create_db_manager('sqlite', db_path='replacements.db')
# For PostgreSQL, use:
# db_manager = create_db_manager('postgres', ini_file='replacements.ini')
mols = list(mutate_mol(m, db_manager, max_size=1))
```

#### Grow Molecule

```python
from rdkit import Chem
from ta_gen.crem.crem import grow_mol
from ta_gen.db import create_db_manager

m = Chem.MolFromSmiles('c1cc(OC)ccc1C')  # toluene
db_manager = create_db_manager('sqlite', db_path='replacements.db')
# For PostgreSQL, use:
# db_manager = create_db_manager('postgres', ini_file='replacements.ini')
mols = list(grow_mol(m, db_manager))
```

#### Link Molecules

```python
from rdkit import Chem
from ta_gen.crem.crem import link_mols
from ta_gen.db import create_db_manager

m1 = Chem.MolFromSmiles('c1cc(OC)ccc1C')  # toluene
m2 = Chem.MolFromSmiles('NCC(=O)O')  # glycine
db_manager = create_db_manager('sqlite', db_path='replacements.db')
# For PostgreSQL, use:
# db_manager = create_db_manager('postgres', ini_file='replacements.ini')
mols = list(link_mols(m1, m2, db_manager))
```

## Examples

See the Jupyter notebooks in the `example/` directory for detailed tutorials:

- `fragmentation_to_db.ipynb`: Demonstrates how to fragment molecules and store results in a database
- `structure_generation.ipynb`: Demonstrates how to generate new molecules using the fragment database

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use oCReM in your research, please cite the paper:

"Expanding accessible chemical space for fragment-based enumeration by orders of magnitude through optimization of the CReM framework"

## Contact

For questions or issues, please contact the authors listed above.