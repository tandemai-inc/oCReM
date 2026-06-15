#!/usr/bin/env python
# -*- coding:utf-8 -*-

import random
import re
from copy import deepcopy
from functools import partial
from multiprocessing import Pool, cpu_count

from rdkit import Chem

from ta_gen.crem_utils.mol_context import combine_core_env_to_rxn_smarts
from ta_gen.utils.fragment_utils import fragment_mol, fragment_mol_link
import multiprocessing
multiprocessing.set_start_method('fork')


def update_protected_ids(mol, protected_ids, replace_ids):
    protected_ids = set(protected_ids) if protected_ids else set()

    if replace_ids:
        ids = set()
        for i in replace_ids:
            ids.update(
                a.GetIdx()
                for a in mol.GetAtomWithIdx(i).GetNeighbors()
                if a.GetAtomicNum() == 1
            )
        ids = (
            set(a.GetIdx() for a in mol.GetAtoms())
            .difference(ids)
            .difference(replace_ids)
        )  # ids which should be protected
        protected_ids.update(
            ids
        )  # since protected_ids has a higher priority add them anyway

    protected_ids = sorted(protected_ids)
    return protected_ids


def zip_new_replacement(new_replacement, input_structure):
    try:
        new_replacement_mol = Chem.MolFromSmiles(new_replacement)
        final_mol = Chem.molzip(input_structure, new_replacement_mol)
        Chem.SanitizeMol(
            final_mol, catchErrors=True
        )  # Bonds can be restored to the aromatic bond type
        # smi = Chem.MolToSmiles(final_mol, isomericSmiles=True)
        smi = Chem.MolToSmiles(Chem.RemoveHs(final_mol), isomericSmiles=True, canonical=True)
        return smi, new_replacement
    except Exception:
        return None, None


def __get_replacements(
    db_manager, env, dist, min_atoms, max_atoms, radius, min_freq=0, **kwargs
):
    condition = [
        f"e.name = '{env}'",
        # f"e.radius = {radius}",
        f"f.core_num_atoms BETWEEN {min_atoms} AND {max_atoms}",
    ]
    if min_freq:
        condition.append(f"ef.frequency >= {min_freq}")
    if isinstance(dist, int):
        condition.append(f"ef.dist2 = {dist}")
    elif isinstance(dist, tuple) and len(dist) == 2:
        condition.append(f"ef.dist2 BETWEEN {dist[0]} AND {dist[1]}")

    sql = f"""
    SELECT f.core_smi
    FROM fragment f
    JOIN env_fragment ef ON f.id = ef.fragment_id
    JOIN env e ON ef.env_id = e.id
    WHERE {' AND '.join(condition)}
    """
    return set(db_manager.execute(sql))


def get_core_smi_replacements(
    db_manager,
    env,
    core_smi,
    dist,
    min_atoms,
    max_atoms,
    max_replacements,
    radius,
    min_freq=0,
    **kwargs,
):
    results = __get_replacements(
        db_manager, env, dist, min_atoms, max_atoms, radius, min_freq, **kwargs
    )
    # remove core_smi from results
    results.discard((core_smi,))

    if max_replacements is not None:
        n = min(len(results), max_replacements)
        selected_results = random.sample(list(results), n)
    else:
        selected_results = list(results)

    for new_core_smi in selected_results:
        yield new_core_smi


def gen_new_replacements(  # noqa: C901
    fragments,
    mol,
    mol_hac,
    db_manager,
    min_inc,
    max_inc,
    max_replacements,
    num_cpus,
    radius,
    min_size=0,
    max_size=8,
    min_rel_size=0,
    max_rel_size=1,
    dist=None,
    min_freq=0,
    products=None,
    return_core=False,
):
    if not products:
        products = set()

    if return_core:
        return_format = lambda _smi, _new_core_smi: (_smi, _new_core_smi)
    else:
        return_format = lambda _smi, _new_core_smi: _smi

    _zip_new_replacement = partial(zip_new_replacement, input_structure=mol)

    with Pool(min(num_cpus, cpu_count())) as p:
        for env_smarts, core_smi, *_ in fragments:
            num_heavy_atoms = Chem.MolFromSmiles(core_smi).GetNumHeavyAtoms()
            hac_ratio = num_heavy_atoms / mol_hac
            min_atoms = num_heavy_atoms + min_inc
            max_atoms = num_heavy_atoms + max_inc

            if not (
                min_size <= num_heavy_atoms <= max_size
                and min_rel_size <= hac_ratio <= max_rel_size
            ):
                continue

            new_replacement_generator = get_core_smi_replacements(
                db_manager,
                env_smarts,
                core_smi,
                dist,
                min_atoms,
                max_atoms,
                max_replacements,
                radius,
                min_freq=min_freq,
            )
            for smi, new_core_smi in p.starmap(
                _zip_new_replacement,
                new_replacement_generator,
                chunksize=100,
            ):
                if smi and smi not in products:
                    products.add(smi)
                    yield return_format(smi, new_core_smi)


def mutate_mol_for_grow(
    mol,
    db_manager,
    radius=3,
    min_atoms=1,
    max_atoms=2,
    dist=None,
    min_freq=0,
    replace_ids=None,
    max_replacements=None,
    protected_ids=None,
    symmetry_fixes=False,
    num_cpus=1,
):
    products = {Chem.MolToSmiles(mol)}

    protected_ids = update_protected_ids(mol, protected_ids, replace_ids)

    f = fragment_mol(
        mol, radius, protected_ids=protected_ids, symmetry_fixes=symmetry_fixes
    )  # [(env smiles, core smiles, list of atom ids)]

    valid_f = []
    for t in f:
        num_heavy_atoms = Chem.MolFromSmiles(t[1]).GetNumHeavyAtoms()
        if num_heavy_atoms == 0:
            valid_f.append(t)

    mol_hac = mol.GetNumHeavyAtoms()
    for env_smarts, core_smi, atom_ids in valid_f:
        for atom_id in atom_ids:
            side_chain = deepcopy(mol)
            side_chain = Chem.AddHs(side_chain)
            atom = side_chain.GetAtomWithIdx(atom_id)
            atom.SetAtomicNum(0)
            atom.SetAtomMapNum(1)
            yield from gen_new_replacements(
                [(env_smarts, core_smi, atom_ids)],
                side_chain,
                mol_hac,
                db_manager,
                min_atoms,
                max_atoms,
                max_replacements,
                num_cpus,
                radius=radius,
                dist=dist,
                min_freq=min_freq,
                min_size=0,
                max_size=0,
                products=products,
            )


def remove_atoms_and_keep_attachments(mol, atoms_to_remove):
    remove_set = set(atoms_to_remove)

    # 1. find bonds to cut
    bonds_to_cut = []
    for bond in mol.GetBonds():
        idx1 = bond.GetBeginAtomIdx()
        idx2 = bond.GetEndAtomIdx()
        if (idx1 in remove_set) != (idx2 in remove_set):
            bonds_to_cut.append(bond.GetIdx())

    # 2. cut bond and add *
    if bonds_to_cut:
        fragmented_mol = Chem.FragmentOnBonds(mol, bonds_to_cut,
                                              dummyLabels=[(0, 0)] * len(bonds_to_cut))
    else:
        fragmented_mol = mol

    # 3. remove atoms in reverse order
    editable_mol = Chem.RWMol(fragmented_mol)
    for idx in sorted(list(remove_set), reverse=True):
        editable_mol.RemoveAtom(idx)

    result_mol = editable_mol.GetMol()
    Chem.SanitizeMol(result_mol)

    # 4. split molecule into frags
    frags = Chem.GetMolFrags(result_mol, asMols=True)

    # 5. filter frags with *
    valid_frags = []
    for frag in frags:
        # check if fragment is dummy only, if not, it is a valid fragment
        is_dummy_only = all(atom.GetAtomicNum() == 0 for atom in frag.GetAtoms())

        if not is_dummy_only:
            valid_frags.append(frag)

    if not valid_frags:
        return ""

    # 6. combine valid frags into a single SMILES string
    combined_smiles = ".".join([Chem.MolToSmiles(f) for f in valid_frags])

    return combined_smiles


def __check_potential_match(env, side_chain, possible_match):
    for env_id, side_chain_id in enumerate(possible_match):
        env_atom = env.GetAtomWithIdx(env_id)
        side_chain_atom = side_chain.GetAtomWithIdx(side_chain_id)
        if env_atom.GetAtomMapNum() > 0 and side_chain_atom.GetAtomicNum() != 0:
            return False

    return True

def __get_potential_macthes(side_chain, env_info):
    matches = set()
    for env_id, _env_info in env_info.items():
        env = _env_info["mol"]
        all_matches = side_chain.GetSubstructMatches(env, uniquify=False)
        for possible_match in all_matches:
            if __check_potential_match(env, side_chain, possible_match):
                matches.add(env_id)
                break

    return matches

def map_env_with_side_chain(env_smarts, side_chain_smi):
    # parse envs
    env_smarts_list = env_smarts.split(".")
    env_info = {}
    for env_smart in env_smarts_list:
        match = re.search(r'\[\*:(\d+)\]', env_smart)
        if match:
            env_id = int(match.group(1))
            env_info[env_id] = {
                "mol": Chem.MolFromSmarts(env_smart), "smart": env_smart
            }


    # parse side chain
    side_chain_smi_list = side_chain_smi.split(".")

    # search potential matched envs for each chain
    side_chain_matched_envs = []
    for side_chain_smi in side_chain_smi_list:
        side_chain = Chem.MolFromSmiles(side_chain_smi)
        matched_envs = __get_potential_macthes(side_chain, env_info)
        side_chain_matched_envs.append([side_chain_smi, matched_envs])

    side_chain_matched_envs = sorted(side_chain_matched_envs, key=lambda x: len(x))
    final_matches = {}
    while side_chain_matched_envs:
        cur_side_chain_info =  side_chain_matched_envs[0]
        side_chain_smi, cur_matched_envs = cur_side_chain_info
        if cur_matched_envs:
            matched_env_id = cur_matched_envs.pop()
            final_matches[side_chain_smi] = matched_env_id
            # update side_chain_matched_envs
            side_chain_matched_envs = side_chain_matched_envs[1:]
            for ele in side_chain_matched_envs:
                ele[1].discard(matched_env_id)
        else:
            final_matches[side_chain_smi] = None
            # update side_chain_matched_envs
            side_chain_matched_envs = side_chain_matched_envs[1:]

    # fill None
    used_env_ids = set([_ for _ in final_matches.values() if _])
    unused_ids = set(range(1, len(env_smarts_list)+1)) - used_env_ids
    for side_chain_smi in final_matches:
        if final_matches[side_chain_smi] is None:
            final_matches[side_chain_smi] = unused_ids.pop()

    # update side chain
    side_chains = []
    for side_chain_smi, env_id in final_matches.items():
        mol = Chem.MolFromSmiles(side_chain_smi)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomMapNum(env_id)
        side_chains.append(Chem.MolToSmiles(mol))

    return ".".join(side_chains)

def mutate_mol(
    mol,
    db_manager,
    radius=3,
    min_inc=-2,
    max_inc=2,
    min_size=0,
    max_size=10,
    min_rel_size=0,
    max_rel_size=1,
    dist=None,
    min_freq=0,
    max_replacements=None,
    replace_ids=None,
    protected_ids=None,
    symmetry_fixes=False,
    num_cpus=1,
):
    products = {Chem.MolToSmiles(mol)}

    protected_ids = update_protected_ids(mol, protected_ids, replace_ids)

    f = fragment_mol(
        mol, radius, protected_ids=protected_ids, symmetry_fixes=symmetry_fixes
    )  # [(env smiles, core smiles, list of atom ids)]

    mol_hac = mol.GetNumHeavyAtoms()
    for env_smarts, core_smi, atom_ids in f:
        if core_smi.count("*") == 1:
            side_chain_smi = remove_atoms_and_keep_attachments(mol, atom_ids)
            side_chain = Chem.MolFromSmiles(side_chain_smi)
            for atom in side_chain.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    atom.SetAtomMapNum(1)
        else:
            side_chain_smi = remove_atoms_and_keep_attachments(mol, atom_ids)
            side_chain_smi = map_env_with_side_chain(env_smarts, side_chain_smi)
            side_chain = Chem.MolFromSmiles(side_chain_smi)

        if not side_chain_smi:
            continue

        yield from gen_new_replacements(
            [(env_smarts, core_smi, atom_ids)],
            side_chain,
            mol_hac,
            db_manager,
            min_inc,
            max_inc,
            max_replacements,
            num_cpus,
            radius=radius,
            dist=dist,
            min_freq=min_freq,
            min_size=min_size,
            max_size=max_size,
            min_rel_size=min_rel_size,
            max_rel_size=max_rel_size,
            products=products,
        )


def grow_mol(
    mol,
    db_manager,
    radius=3,
    min_atoms=1,
    max_atoms=2,
    max_replacements=None,
    replace_ids=None,
    symmetry_fixes=False,
    num_cpus=1,
):
    mol = Chem.AddHs(mol)
    protected_ids = set()
    if replace_ids:
        ids = set()  # ids if replaceable Hs
        for i in replace_ids:
            if mol.GetAtomWithIdx(i).GetAtomicNum() == 1:
                ids.add(i)
            else:
                ids.update(
                    a.GetIdx()
                    for a in mol.GetAtomWithIdx(i).GetNeighbors()
                    if a.GetAtomicNum() == 1
                )
        ids = set(
            a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1
        ).difference(
            ids
        )  # ids of Hs to protect
        protected_ids.update(ids)  # Hs should be protected

    return mutate_mol_for_grow(
        mol,
        db_manager,
        radius,
        min_atoms=min_atoms,
        max_atoms=max_atoms,
        max_replacements=max_replacements,
        replace_ids=None,
        protected_ids=protected_ids,
        num_cpus=num_cpus,
        symmetry_fixes=symmetry_fixes,
    )


def mark_wildcard_by_env(mol, env):
    matches = mol.GetSubstructMatches(env, uniquify=False)
    if not matches:
        return False

    wildcard_indices_in_env = [
        a.GetIdx()
        for a in env.GetAtoms()
        if (a.GetAtomicNum() == 0) & (a.GetAtomMapNum() != 0)
    ]

    for match in matches:
        for env_idx in wildcard_indices_in_env:
            mol_atom_idx = match[env_idx]
            atom = mol.GetAtomWithIdx(mol_atom_idx)
            if atom.GetAtomicNum() == 0:
                return True
    return False


def combine_link_mols(side_chain_1, side_chain_2, env_smarts):
    mol1 = Chem.MolFromSmiles(Chem.MolToSmiles(side_chain_1))
    mol2 = Chem.MolFromSmiles(Chem.MolToSmiles(side_chain_2))

    env1, env2 = env_smarts.split(".")
    for env in env_smarts.split("."):
        if "*:1" in env:
            env1 = env
        elif "*:2" in env:
            env2 = env

    env1 = Chem.MolFromSmarts(env1)
    env2 = Chem.MolFromSmarts(env2)

    mol1_in_env_1 = mark_wildcard_by_env(mol1, env1)
    mol2_in_env_2 = mark_wildcard_by_env(mol2, env2)
    if mol1_in_env_1 and mol2_in_env_2:
        for atom in mol2.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomMapNum(2)

        combined = Chem.CombineMols(mol1, mol2)
    else:
        for atom in mol1.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomMapNum(2)
        for atom in mol2.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomMapNum(1)
        combined = Chem.CombineMols(mol1, mol2)
    return combined


def link_mols(
    mol1,
    mol2,
    db_manager,
    radius=3,
    min_atoms=1,
    max_atoms=2,
    dist=None,
    min_freq=0,
    max_replacements=None,
    replace_ids_1=None,
    replace_ids_2=None,
    num_cpus=1,
):
    mol1 = Chem.AddHs(mol1)
    mol2 = Chem.AddHs(mol2)

    protected_ids_1 = set()
    if replace_ids_1:
        ids = set()  # ids if replaceable Hs
        for i in replace_ids_1:
            if mol1.GetAtomWithIdx(i).GetAtomicNum() == 1:
                ids.add(i)
            else:
                ids.update(
                    a.GetIdx()
                    for a in mol1.GetAtomWithIdx(i).GetNeighbors()
                    if a.GetAtomicNum() == 1
                )
        ids = set(
            a.GetIdx() for a in mol1.GetAtoms() if a.GetAtomicNum() == 1
        ).difference(
            ids
        )  # ids of Hs to protect
        protected_ids_1.update(ids)  # Hs should be protected

    protected_ids_2 = set()
    if replace_ids_2:
        ids = set()  # ids if replaceable Hs
        for i in replace_ids_2:
            if mol2.GetAtomWithIdx(i).GetAtomicNum() == 1:
                ids.add(i)
            else:
                ids.update(
                    a.GetIdx()
                    for a in mol2.GetAtomWithIdx(i).GetNeighbors()
                    if a.GetAtomicNum() == 1
                )
        ids = set(
            a.GetIdx() for a in mol2.GetAtoms() if a.GetAtomicNum() == 1
        ).difference(
            ids
        )  # ids of Hs to protect
        protected_ids_2.update(ids)  # Hs should be protected

    fragments = fragment_mol_link(
        mol1,
        mol2,
        radius,
        protected_ids_1=protected_ids_1,
        protected_ids_2=protected_ids_2,
    )  # [(env smiles, core smiles, list of atom ids)]

    for env_smarts, core_smi, atom_ids_1, atom_ids_2 in fragments:
        side_chain_smi_1 = remove_atoms_and_keep_attachments(mol1, atom_ids_1)
        # replaced_string = re.sub(r"\*", r"\[*:1\]", side_chain_smi_1)  # [1*] -> [*:1]
        side_chain_1 = Chem.MolFromSmiles(side_chain_smi_1)
        for atom in side_chain_1.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomMapNum(1)

        side_chain_smi_2 = remove_atoms_and_keep_attachments(mol2, atom_ids_2)
        # replaced_string = re.sub(r"\*", r"\[*:1\]", side_chain_smi_2)  # [1*] -> [*:1]
        side_chain_2 = Chem.MolFromSmiles(side_chain_smi_2)
        for atom in side_chain_2.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomMapNum(2)

        mol = combine_link_mols(side_chain_1, side_chain_2, env_smarts)
        mol_hac = mol.GetNumHeavyAtoms()

        yield from gen_new_replacements(
            [(env_smarts, core_smi, atom_ids_1, atom_ids_2)],
            mol,
            mol_hac,
            db_manager,
            min_atoms,
            max_atoms,
            max_replacements,
            num_cpus,
            radius,
            dist=dist,
            min_freq=min_freq,
            min_size=0,
            max_size=0,
        )


if __name__ == "__main__":
    x = map_env_with_side_chain(
        "*-C-O-C(=O)-[*:1].*1:c:c:c(-[*:3]):c:c:1.O-[*:2]",
        "*C(=O)OCC.*O.*c1ccccc1"
    )
    print(x)
