#!/bin/python

import re
from multiprocessing import Pool
from typing import List

from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition

from ta_gen.utils.logger import LOGGER


def update_attachment_point_indexes(scaffolds):
    """Update attachment point indexes - make it start from 1
    Args:
        scaffolds (list): annotated smiles list
    Returns:
        list: updated annotated smiles list
    """
    pattern = r"\[\*\:(\d+)\]"
    replacement = "[*:{new_idx}]"
    updated_scaffolds = []
    for s in scaffolds:
        # match.group(0) -> whole matched string -> [*:0]
        # match.group(1) -> first sub group -> (\d+) -> 0
        new_s = re.sub(
            pattern, lambda match: replacement.format(new_idx=int(match.group(1)) + 1), s
        )
        new_s = new_s.replace("\\", "").replace("/", "")  # remove directionality of chemical bonds

        mol = Chem.MolFromSmiles(new_s)
        rm = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                rm.extend([a.GetIdx() for a in atom.GetNeighbors()])

        for i in rm:
            atom = mol.GetAtomWithIdx(i)
            atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)

        smiles_no_chirality = Chem.MolToSmiles(mol)

        Chem.SanitizeMol(mol)
        updated_scaffolds.append(smiles_no_chirality)
    return updated_scaffolds


def run_rgroup_decompose(smi: str, scaffolds: List[str]):
    core_mols = [Chem.MolFromSmiles(s) for s in scaffolds]
    ps = Chem.AdjustQueryParameters.NoAdjustments()
    ps.makeDummiesQueries = True

    query_core_mols = [Chem.AdjustQueryProperties(c, ps) for c in core_mols]

    mol = Chem.MolFromSmiles(smi.replace("\\", "").replace("/", ""))
    Chem.AddHs(mol, addCoords=True)

    groups, _ = rdRGroupDecomposition.RGroupDecompose(
        query_core_mols, [mol], asSmiles=True, asRows=True
    )
    return groups


def match_mols_with_cores(smi_list, mols, query_core_mols, core_names):
    matched_mols = []
    core_name_list = []
    unmatched_smis = []
    for smi, m in zip(smi_list, mols):
        find = False
        for i, qc in enumerate(query_core_mols):
            if m.HasSubstructMatch(qc):
                matched_mols.append(m)
                core_name_list.append(core_names[i])
                find = True
                break
        if not find:
            unmatched_smis.append(smi)
            print(f"WARNING: No matching core found for SMILES: {smi}")
    return matched_mols, core_name_list, unmatched_smis


def identify_rgroups(
    smi_list,
    scaffolds,
    attachment_points_mapping,
    num_cpus=1,
):
    """Process of R-Groups identification
    Args:
        smi_list (list): generated smiles
        scaffolds (list): annotated smiles list
        attachment_points_mapping (dict): mapping of indexes
        chunk_size (int): chunk size of rdRGroupDecomposition.RGroupDecompose
        num_cpus (int): number of cpus for multiprocessing

    Returns:
        list: identified R-Groups
    """
    core_names = []
    core_name_scaffold_dict = dict()
    attachment_point_cnt = dict()
    for i, smi in enumerate(scaffolds):  # scaffold start from 0
        name = f"core_{i}"
        core_names.append(name)
        core_name_scaffold_dict[name] = smi
        attachment_point_cnt[name] = smi.count("*")

    scaffolds = update_attachment_point_indexes(scaffolds)  # scaffold start from 1

    # convert smiles string to mol object and store in a list
    mols = [Chem.MolFromSmiles(smi.replace("\\", "").replace("/", "")) for smi in smi_list]
    # add hydrogens to all molecules
    mols = [Chem.AddHs(m, addCoords=True) for m in mols]

    # create list of mol object from core SMILES
    core_mols = [Chem.MolFromSmiles(s) for s in scaffolds]  # core rdmol list

    # set atom properties to use for next steps
    # Returns an AdjustQueryParameters object with all parameters set to false
    ps = Chem.AdjustQueryParameters.NoAdjustments()
    ps.makeDummiesQueries = (
        True  # convert dummy atoms(*) without isotope labels to any-atom queries
    )

    # make a list of all query cores
    # Chem.AdjustQueryProperties Returns a new molecule
    # where the query properties of atoms have been modified
    query_core_mols = [Chem.AdjustQueryProperties(c, ps) for c in core_mols]

    # perform substructure match and if matches, collect it in a new list named 'mms'
    matched_mols, core_name_list, unmatched_smis = match_mols_with_cores(
        smi_list, mols, query_core_mols, core_names
    )

    if len(mols) != len(matched_mols):
        print("WARNING: Not all Substructure matches are found")
        print("Number of unmatched mols:", len(unmatched_smis))
        smi_list = [s for s in smi_list if s not in unmatched_smis]

    # an empty list that will collect r-group decomposition values for all molecules
    all_groups = []

    inputs = []
    for smi in smi_list:
        inputs.append((smi, scaffolds))
    LOGGER.info(f"run RGroupDecompose on {num_cpus} cpus... ")
    with Pool(num_cpus) as pool:
        res_list = pool.starmap(run_rgroup_decompose, inputs)
    for res in res_list:
        all_groups.extend(res)

    # [{'Core': 'O=C1NC2(CC([*:1])C2)C(=O)N1[*:2]', 'R1': 'ClCC[*:1]', 'R2': '[H][*:2]'},
    #  {'Core': 'O=C1NC2(CC([*:1])C2)C(=O)N1[*:2]', 'R1': '[H][*:1]', 'R2': 'C=C(Br)C[*:2]'}]
    valid_smi_list = []
    rgroups = []
    for i, res in enumerate(all_groups):
        if not res:
            continue
        valid_smi_list.append(smi_list[i])
        core_name = core_name_list[i]
        rg = [""] * attachment_point_cnt[core_name]
        for k, v in res.items():
            # only extract dictionary items that have heavy atom as R group
            if k == "Core" or "[H]" in v:
                continue
            idx = int(k[1:]) - 1  # R2 -> 1
            if 0 <= idx < len(rg):
                # rg[idx] = v
                # convert v to original attachment point
                # R2: C=C(Br)C[*:2]
                # 2 -> 1 -> new2old {1: 0} -> 0
                # C=C(Br)C[*:0]
                scaffold = core_name_scaffold_dict[core_name]
                old_idx = attachment_points_mapping[scaffold]["new2old"][idx]
                rg[idx] = v.replace(f"[*:{idx + 1}]", f"[*:{old_idx}]")

        rgroups.append(";".join(rg))

    LOGGER.info("RGroupDecompose done")
    return smi_list, rgroups
