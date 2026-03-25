# oCReM

This repository provides the database and tools for the paper "Expanding accessible chemical space for fragment-based enumeration by orders of magnitude through optimization of the CReM framework".

It contains the data of modified CrEM library and the original ChEMBL fragment library, along with tools for molecular fragmentation, fragment database management, and structure generation.

## Authors

- Shaojin Hu
- Qinyu Chen
- Yinhui Yi
- Paul Pilot
- James Xu
- Abir Ganguly
- Albert C. Pan

## Data Access

- **Original ChemBL fragment library**: [http://www.qsar4u.com/pages/crem.php](https://www.qsar4u.com/pages/crem.php)
- **Download the oCReM databases**:
  - **ChEMBL22**:
    - [SQLite](https://zenodo.org/records/19107922/files/ocrem_chembl22_sqlite.tar.gz)
    - [PostgreSQL](https://zenodo.org/records/19107922/files/ocrem_chembl22_postgres.dump)
  - **ChEMBL36**:
    - [SQLite](https://zenodo.org/records/19107922/files/ocrem_chembl_36_sqlite.tar.gz)
    - [PostgreSQL](https://zenodo.org/records/19107922/files/ocrem_chembl36_postgres.dump)
## Overview

oCReM (optimized Chemical Reassembly) is a framework for fragment-based molecular design that includes:

- **Molecular Fragmentation**: Breaking down molecules into structural fragments
- **Fragment Database Management**: Storing and retrieving fragments in SQLite or PostgreSQL databases
- **Structure Generation**: Generating new molecules through fragment-based assembly
  - **Mutate**: Replace fragments in a molecule
  - **Grow**: Extend a molecule with new fragments
  - **Link**: Connect multiple molecules using linker fragments

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

### Clone the Repository

First, clone the oCReM repository:

```bash
git clone https://github.com/your-username/oCReM.git
cd oCReM
```

### Create Conda Environment

Make sure you have installed Miniconda.

```bash
conda env create -f environment.yml
conda activate crem
# Set Python path to include the project root
export PYTHONPATH=$PYTHONPATH:$(pwd)
```

### 1. Molecular Fragmentation

#### Input File Format

The input file (`test.smi`) should contain one molecule per line in SMILES format:

```
CC(=O)Oc1ccccc1C(=O)O
CCO
c1ccccc1
```

#### Fragment to File
Generate fragments for each molecule in the input file and save them to a CSV file.

```bash
python ta_gen/bin/fragmentation.py --input test.smi --out test_frag.csv --mode 0 --ncpu 10 --radius 3
```

#### Fragment to SQLite Database
Generate fragments for each molecule in the input file and save them to a SQLite database. If the SQLite database file does not exist, it will be automatically created.

```bash
python ta_gen/bin/fragmentation.py --input test.smi --mode 0 --ncpu 10 --radius 3 --use_db --db_type sqlite --db_path test.db
```

#### Fragment to PostgreSQL Database

##### PostgreSQL Configuration File

Create a configuration file (`test.ini`) with the following format:

```ini
[database]
host=your_host
port=your_port
user=your_username
password=your_password
database=your_database
```

##### Command
Generate fragments for each molecule in the input file and save them to a PostgreSQL database. If the database does not exist, it will try to create it.

```bash
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

### 3. Using the oCReM Databases

#### SQLite Database

1. Download the SQLite database tar.gz file (e.g., `ocrem_chembl36_sqlite.tar.gz`)
2. Extract the tar.gz file:

```bash
tar -xzf ocrem_chembl36_sqlite.tar.gz
```

This will extract to a file named `chembl_36.db`.

3. Use the extracted SQLite database file with the `create_db_manager` function:

```python
from ta_gen.db import create_db_manager

db_manager = create_db_manager('sqlite', db_path='path/to/chembl_36.db')
```

#### PostgreSQL Database

**Note: Make sure you have PostgreSQL installed before proceeding.**

1. Download the PostgreSQL dump file (e.g., `ocrem_chembl36_postgres.dump`)
2. Create a new PostgreSQL database:

```bash
createdb -U your_username your_database
# With host and port specified:
# createdb -h your_host -p your_port -U your_username your_database
```

3. Import the dump file into your database:

```bash
pg_restore -U your_username -d your_database ocrem_chembl36_postgres.dump
# With host and port specified:
# pg_restore -h your_host -p your_port -U your_username -d your_database ocrem_chembl36_postgres.dump
```

4. Create PostgreSQL database configuration file:

Create a `database.ini` file with the following format:

```ini
[database]
host=your_host
port=your_port
user=your_username
password=your_password
database=your_database
```

5. Use it with the `create_db_manager` function:

```python
from ta_gen.db import create_db_manager

db_manager = create_db_manager('postgres', ini_file='path/to/database.ini')
```

## Examples

See the Jupyter notebooks in the `example/` directory for detailed tutorials:

- **fragmentation_to_db.ipynb**: Demonstrates how to fragment molecules and store results in a database
- **structure_generation.ipynb**: Demonstrates how to generate new molecules using the fragment database

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use oCReM in your research, please cite the paper:

> "Expanding accessible chemical space for fragment-based enumeration by orders of magnitude through optimization of the CReM framework"

## Contact

For questions or issues, please contact the authors listed above.