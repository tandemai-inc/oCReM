#!/bin/python
import os
import re
from collections import Counter
from pathlib import Path
from typing import List

import pandas as pd
from rdkit import Chem

from ta_gen.utils.common_utils import standardize_smiles

ALERTS = os.path.join(Path(__file__).parent.absolute(), "alerts.csv")
rule_df = pd.read_csv(ALERTS)
df = rule_df.dropna()
rules_list = []
for _, row in df.iterrows():
    smarts = row["smarts"]
    pattern = Chem.MolFromSmarts(smarts)
    rule_id = row["rule_id"]
    if pattern:
        rules_list.append((rule_id, pattern, row["max"]))
    else:
        print(f"Error parsing SMARTS for rule {rule_id}")


def generate_tafilter_id(smiles: str, scaffold: str = None) -> List:
    """Generate tafilter id

    Args:
        smiles (str): molecule smiles.

    Returns:
        list: matched rule id list.

    """
    tafilter_id = []
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return tafilter_id

    if scaffold:
        # remove [*:d], [*]
        scaffold = re.sub(r"\[\*:(\d+)\]", "", scaffold).replace("[*]", "").replace("()", "")
        scaffold_smart = Chem.MolFromSmarts(standardize_smiles(scaffold))
        mapping = mol.GetSubstructMatches(scaffold_smart)
        scaffold_matches = set()
        for t in mapping:
            for i in t:
                scaffold_matches.add(i)
            break  # if matches more than one time, only keep one

    for rule_id, pattern, max_val in rules_list:
        structure_matches = mol.GetSubstructMatches(pattern)
        if scaffold:
            structure_matches = [t for t in structure_matches if not set(t) <= scaffold_matches]
        if len(structure_matches) > max_val:
            tafilter_id.append(rule_id)
    return tafilter_id


def get_scaffold_pattern_count(scaffold: str = None):
    cnt = Counter()
    if scaffold:
        scaffold = re.sub(r"\[\*:(\d+)\]", "[H]", scaffold).replace("[*]", "[H]").replace("()", "")
        scaffold_mol = Chem.MolFromSmiles(standardize_smiles(scaffold))
        for rule_id, pattern, _ in rules_list:
            match_cnt = len(scaffold_mol.GetSubstructMatches(pattern))
            cnt[rule_id] = match_cnt
    return cnt


def get_side_chain_pattern_count(side_chain_smi: str = None):
    cnt = Counter()
    if side_chain_smi:
        side_chain_mol = Chem.MolFromSmiles(side_chain_smi)
        for rule_id, pattern, _ in rules_list:
            match_cnt = len(side_chain_mol.GetSubstructMatches(pattern))
            cnt[rule_id] = match_cnt
    return cnt


def generate_tafilter_id_by_count(
    smiles: str,
    cnt_scaffold_pattern: dict,
    cnt_side_chain_pattern: dict,
    ignore_tf_alerts: set = None,
) -> List:
    """Generate tafilter id
    Args:
        smiles (str): molecule smiles.
        cnt_scaffold_pattern (dict): scaffold pattern count.
        side_chain_smi (dict): side chain pattern count.
        ignore_tf_alerts (set): ignore rule id set.

    Returns:
        list: matched rule id list.
    """
    ignore_tf_alerts = ignore_tf_alerts or set()

    tafilter_id = []
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return tafilter_id

    for rule_id, pattern, max_val in rules_list:
        if rule_id in ignore_tf_alerts:
            continue
        structure_matches = mol.GetSubstructMatches(pattern)
        count = (
            len(structure_matches)
            - cnt_scaffold_pattern.get(rule_id, 0)
            - cnt_side_chain_pattern.get(rule_id, 0)
        )
        if count > max_val:
            tafilter_id.append(rule_id)
    return tafilter_id


def get_sdf_tafilter_id(sdf: str) -> set:
    mol = Chem.MolFromMolFile(sdf)
    if not mol:
        return set()
    tafilter_id = set()
    for rule_id, pattern, max_val in rules_list:
        structure_matches = mol.GetSubstructMatches(pattern)
        if len(structure_matches) > max_val:
            tafilter_id.add(rule_id)
    return tafilter_id
