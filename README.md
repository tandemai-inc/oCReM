# oCrEM

This repository provides the database and tools for the paper "Expanding accessible chemical space for fragment-based enumeration by orders of magnitude through optimization of the CReM framework". 

It contains the data of modified CrEM library and the original ChEMBL fragment library.

## Authorlist

- Shaojin Hu, Qinyu Chen, Yinhui Yi, Paul Pilot, James Xu, Abir Ganguly and Albert C. Pan

## Data Access
Original ChemBL fragment library: [http://www.qsar4u.com/pages/crem.php](https://www.qsar4u.com/pages/crem.php)

## Create conda environment
Make sure you have installed miniconda.
```bash
conda env create -f environment.yml
```

## Create Postgres Database
Make sure you have installed postgres and created a database.

### Download the modified CrEM library
You can download the modified CrEM library from [https://zenodo.org/records/17796862/files/crem.tar.gz?download=1](https://zenodo.org/records/17796862/files/crem.tar.gz?download=1).
```
mkdir run
cd run
wget https://zenodo.org/records/17796862/files/crem.tar.gz
```

### Unzip the modified CrEM library
```bash
tar -xzf crem.tar.gz
```

### Import the modified CrEM library to the postgres database
You should create .ini file for postgres database, the template is ./import_db/db.ini
```bash
cp ../import_db/db.ini db.ini
# Edit the db.ini file to set the database connection parameters
# activate conda env
conda activate crem
# Run the import script
python ../import_db/import_db_postgres.py -i db.ini -d ./data
```

## Usage

### Create config file

The template is under config_template/config.yaml. Set the parameter-crem-db_config to db.ini.

### Run
```bash
# copy config.yaml to run directory
cp ../config_template/config.yaml config.yaml
# create scaffolds_file and set smart in it, for example:
# scaffolds_file: input.smi
# the content in scaffolds_file:
# [*:1]c1ccnc(N[*:2])c1
../ta_gen/bin/TaGEN -i config.yaml
```



