#!/usr/bin/env python
# -*- coding:utf-8 -*-

import ast
import configparser
import functools
import inspect
import re
import traceback
from pathlib import Path
from typing import Iterator, List

import pandas as pd
from rdkit import Chem
from rdkit.Chem import QED, AllChem, Descriptors, rdFreeSASA, rdMolDescriptors

from ta_gen.utils.const import MAXINUM_NUM_OF_OUTPUT_MOLS

from .sascore import calculateScore


def dict_to_cmdline(argdict):
    """A utility function that transforms an arg=value dictionary to a command-line option string.

    Args:
        argdict (dict): Dictionary whose key-value pairs represent command-line
            argument names and values.

    Returns:
        A string that reproduces the string one will type when running the
            script.
    """
    cmd = []
    for key, value in argdict.items():
        arg = "--" + str(key)
        if type(value) is bool:
            if value is True:
                cmd.append(arg)
        elif type(value) in [list, set, tuple]:
            cmd.append(arg)
            for v in value:
                cmd.append(str(v))
        elif type(value) in [dict]:
            cmd.append(arg)
            for k, v in value.items():
                cmd.extend([str(k), str(v)])
        else:
            cmd.append(arg)
            cmd.append('"' + str(value).replace('"', '\\"') + '"')

    return " ".join(cmd)


def standardize_smiles(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def validate_smiles(smiles: str) -> bool:
    """
    Validates a SMILES string using RDKit.

    Args:
        smiles (str): A SMILES string to validate.

    Returns:
        bool: True if the SMILES string is valid, False otherwise.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        else:
            return True
    except Exception:
        return False


def get_db_config(config_path="config.ini"):
    config = configparser.ConfigParser()
    config.read(config_path)
    db_config = {
        "host": config["database"]["host"],
        "port": config["database"]["port"],
        "user": config["database"]["user"],
        "password": config["database"]["password"],
        "dbname": config["database"]["database"],
    }
    return db_config


def remove_mark(smiles: str) -> str:
    """
    Removes the *:n mark from a SMILES string and replaces it with an H atom.

    Args:
        smiles (str): The SMILES string to modify.

    Returns:
        str: The modified SMILES string.
    """
    regex = re.compile(r"\*:\d+")
    return regex.sub("H", smiles)


def validate_scaffold(scaffold: str):
    """
    Validates a scaffold string using RDKit.

    Args:
        scaffold (str): A scaffold string to validate.

    Returns:
        bool: True if the scaffold string is valid, False otherwise.
    """
    pattern = r"\[\*\:(\d+)\]"
    matches = re.findall(pattern, scaffold)
    if not matches:  # no attachment point
        return False
    if len(matches) != len(set(matches)):  # duplicated attachment points
        return False
    smiles = remove_mark(scaffold)
    return validate_smiles(smiles)


def calculate_num_ionizable_groups(mol, ionizable_structure_smarts) -> int:
    """
    Calculates the number of ionizable groups in a compound.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object.
        ionizable_structure_smarts (pd.Series): List of ionizable structure SMARTS strings.

    Returns:
        int: Number of ionizable groups.
    """
    try:
        substructure_matches = ionizable_structure_smarts.apply(
            lambda x: mol.HasSubstructMatch(Chem.MolFromSmarts(x))
        )
        num_ionizable_groups = sum(substructure_matches)
    except Exception as e:
        print(f"Error in calculate_num_ionizable_groups: {e}")
        num_ionizable_groups = -1
    return num_ionizable_groups


def calculate_mol_properties(smiles: str, ionizable_structure_smarts) -> dict:
    """
    Calculates molecular properties of a compound from its SMILES string.

    Args:
        smiles (str): SMILES string of the compound
        ionizable_structure_smarts (pd.Series): List of ionizable structure SMARTS strings.

    Returns:
        dict: molecular weight (MW), polar surface area (PSA), hydrogen bond acceptors (HBA),
        hydrogen bond donors (HBD), LogP, rotatable bonds (rotB), number of atoms (numA),
        and molar refractivity (mr).
    """
    # Convert SMILES string to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {
            "MW": -999,
            "PSA": -999,
            "HBA": -999,
            "HBD": -999,
            "LogP": -999,
            "rotB": -999,
            "numA": -999,
            "mr": -999,
            "sascore": -999,
            "numAromaticRings": 0,
            "formula": "",
            "qed": -999,
            "numChiralCenters": 0,
            "numIonizableGroups": 0,
        }

    # Calculate properties using RDKit descriptors
    mw = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)
    hba = Descriptors.NumHAcceptors(mol)
    hbd = Descriptors.NumHDonors(mol)
    clogp = Descriptors.MolLogP(mol)
    rotB = Descriptors.NumRotatableBonds(mol)
    numA = mol.GetNumHeavyAtoms()
    mr = Descriptors.MolMR(mol)
    sascore = calculateScore(mol)
    num_aromatic_rings = Descriptors.NumAromaticRings(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    qed = QED.qed(mol)
    num_chiral_centers = len(
        Chem.FindMolChiralCenters(
            mol, force=True, includeUnassigned=True, includeCIP=True
        )
    )
    num_ionizable_groups = calculate_num_ionizable_groups(
        mol, ionizable_structure_smarts
    )

    # Return dictionary of calculated properties
    return {
        "MW": round(mw, 4),
        "PSA": round(tpsa, 4),
        "HBA": int(hba),
        "HBD": int(hbd),
        "LogP": round(clogp, 4),
        "rotB": int(rotB),
        "numA": int(numA),
        "mr": round(mr, 4),
        "sascore": round(sascore, 4),
        "numAromaticRings": int(num_aromatic_rings),
        "formula": formula,
        "qed": round(qed, 4),
        "numChiralCenters": int(num_chiral_centers),
        "numIonizableGroups": num_ionizable_groups,
    }


def get_random_sample_from_df(df, num=MAXINUM_NUM_OF_OUTPUT_MOLS):
    if num == -1:
        return df
    total_records = len(df)
    if total_records > num:
        print(
            f"Num of rows {total_records} > {num}, will return a random sample({num}) of items"
        )
        df = df.sample(n=num)
    return df


def load_ionizable_structure_smarts():
    """
    Loads the ionizable structure SMARTS file.

    Returns:
        list: List of ionizable structure SMARTS strings.
    """
    df = pd.read_csv(
        Path(__file__).parents[1] / "data/smarts_rulebook_ph_range.txt", sep=r"\s+"
    )
    return df["SMARTS_TSI"]


def get_smiles_list_from_csv(csv_file):
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"{csv_file} not found")
    df = pd.read_csv(csv_file)
    if "smiles" not in df.columns:
        raise ValueError(f"smiles column not found in {csv_file}")
    smi_list = df["smiles"].to_list()
    return smi_list


def load_db_query(db_query):
    if db_query:
        db_query = ast.literal_eval(db_query)
    else:
        db_query = {}
    return db_query


def read_scaffolds(scaffolds_file: str) -> List[str]:
    """
    Reads a file of scaffolds and returns a list of the scaffold strings.

    Args:
        scaffolds_file (str): The path to the file containing the scaffolds.

    Returns:
        List[str]: A list of scaffold strings.
    """
    scaffolds = []
    with open(scaffolds_file, "r") as f:
        input_lines = f.readlines()
    for s in input_lines:
        s = s.strip()
        if s.startswith("#"):
            continue
        scaffold = s.strip()
        scaffolds.append(scaffold)
    return scaffolds


def fault_tolerance(func):

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception:
            sig = inspect.signature(func)
            traceback_msg = traceback.format_exc()
            try:
                bound_args = sig.bind(*args, **kwargs)
                bound_args.apply_defaults()  # apply default values
                params_str = ", ".join(
                    [
                        f"{name}={repr(value)}"
                        for name, value in bound_args.arguments.items()
                    ]
                )

                print(
                    f"Failed run {func.__name__} with input ({params_str})\n{traceback_msg}"
                )
            except Exception:
                print(f"Failed run {func.__name__} and get inputs.\n{traceback_msg}")

    return wrapper
