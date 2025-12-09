# CrEM-Database

This repository provides the database and tools for the paper "[Your Paper Title Here]". It contains the data of modified CrEM library and the original ChEMBL fragment library.


## Data Access
Original ChemBL fragment library: [http://www.qsar4u.com/pages/crem.php](https://www.qsar4u.com/pages/crem.php)


## Create Postgres Database

### Download the modified CrEM library
You can download the modified CrEM library from [https://zenodo.org/records/17796862/files/crem.tar.gz?download=1](https://zenodo.org/records/17796862/files/crem.tar.gz?download=1).

### Unzip the modified CrEM library
```bash
tar -xzf crem.tar.gz
```

### Import the modified CrEM library to the postgres database
You should create .ini file for postgres database, the template is ./import_db/db.ini
```bash
python import_db/import_db_postgre.py -i frag.ini -d ./data
```

## Create conda environment
```bash
conda env create -f environment.yml
```

## Usage

### Create config file

The template is under ta_gen/config.yml. Set the db_config to frag.ini.

### Run
```bash
./ta_gen/bin/TaGEN -i config.yml
```



