#!/usr/bin/env python
# -*- coding:utf-8 -*-

import re
from copy import deepcopy
from rdkit import Chem
from multiprocessing import Pool, cpu_count
import random
from functools import partial


from ta_gen.utils.fragment_utils import fragment_mol, fragment_mol_link
from ta_gen.crem_utils.mol_context import combine_core_env_to_rxn_smarts

def update_protected_ids(mol, protected_ids, replace_ids):
    protected_ids = set(protected_ids) if protected_ids else set()

    if replace_ids:
        ids = set()
        for i in replace_ids:
            ids.update(
                a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors() if a.GetAtomicNum() == 1
            )
        ids = (
            set(a.GetIdx() for a in mol.GetAtoms()).difference(ids).difference(replace_ids)
        )  # ids which should be protected
        protected_ids.update(ids)  # since protected_ids has a higher priority add them anyway

    protected_ids = sorted(protected_ids)
    return protected_ids


def zip_new_replacement(new_replacement, input_structure):
    try:
        new_replacement_mol = Chem.MolFromSmiles(new_replacement)
        final_mol = Chem.molzip(input_structure, new_replacement_mol)
        Chem.SanitizeMol(
            final_mol, catchErrors=True
        )  # Bonds can be restored to the aromatic bond type
        smi = Chem.MolToSmiles(final_mol, isomericSmiles=True)
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        return smi, new_replacement
    except Exception:
        return None, None


def __get_replacements(
    db_manager, env, dist, min_atoms, max_atoms, radius, min_freq=0, **kwargs
):
    condition = [
        f"e.name = '{env}'",
        f"e.radius = {radius}",
        f"ef.frequency >= {min_freq}",
        f"f.core_num_atoms BETWEEN {min_atoms} AND {max_atoms}",
    ]
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
    db_manager, env, core_smi, dist, min_inc, max_inc, max_replacements, radius, min_freq=0, **kwargs
):
    num_heavy_atoms = Chem.MolFromSmiles(core_smi).GetNumHeavyAtoms()
    min_atoms = num_heavy_atoms + min_inc
    max_atoms = num_heavy_atoms + max_inc

    results = __get_replacements(
        db_manager, env, dist, min_atoms, max_atoms, radius, min_freq, **kwargs
    )

    if max_replacements is not None:
        n = min(len(results), max_replacements)
        selected_results = random.sample(list(results), n)
    else:
        selected_results = list(results)

    for new_core_smi in selected_results:
        if new_core_smi != core_smi:
            yield new_core_smi


def gen_new_replacements(  # noqa: C901
    fragments,
    mol,
    db_manager,
    min_inc,
    max_inc,
    max_replacements,
    num_cpus,
    radius,
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
            new_replacement_generator = get_core_smi_replacements(
                db_manager,
                env_smarts,
                core_smi,
                dist,
                min_inc,
                max_inc,
                max_replacements,
                radius,
                min_freq=min_freq,
            )
            for smi, new_core_smi in p.starmap(
                _zip_new_replacement, new_replacement_generator, chunksize=100,
            ):
                if smi and smi not in products:
                    products.add(smi)
                    yield return_format(smi, new_core_smi)


def mutate_mol(
    mol,
    db_manager,
    radius=3,
    min_inc=-2,
    max_inc=2,
    dist= None,
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

    print("valid f", valid_f)
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
                db_manager,
                min_inc,
                max_inc,
                max_replacements,
                num_cpus,
                radius=radius,
                dist=dist,
                min_freq=min_freq,
                products=products,
            )


def mutate_mol_ch(
    mol,
    db_manager,
    radius=3,
    min_inc=-2,
    max_inc=2,
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

    print(f)

    for env_smarts, core_smi, atom_ids in f:
        print(env_smarts, core_smi, atom_ids)

        smarts = combine_core_env_to_rxn_smarts(core_smi, env_smarts, keep_h=False)

        smarts_mol = Chem.MolFromSmarts(smarts)

        total_match = mol.GetSubstructMatch(smarts_mol)

        core_mol = Chem.MolFromSmarts(core_smi)

        dummy_idx_bond_idx_map = {}
        for atom in core_mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                neighbor = atom.GetNeighbors()
                dummy_idx_bond_idx_map[atom.GetAtomMapNum()] = [
                    atom.GetIdx(),
                    neighbor[0].GetIdx(),
                ]

        core_matches = mol.GetSubstructMatches(core_mol)

        frag_mol = None

        for core_match in core_matches:
            if set(core_match).issubset(set(total_match)):
                cut_bonds = []
                labels = []
                for k, v in dummy_idx_bond_idx_map.items():
                    cut_bonds.append(
                        mol.GetBondBetweenAtoms(core_match[v[0]], core_match[v[1]]).GetIdx()
                    )
                    labels.append((k, k))
                frag_mol = Chem.FragmentOnBonds(
                    mol, cut_bonds, addDummies=True, dummyLabels=labels
                )
                break

        if not frag_mol:
            continue

        frag_smi = Chem.MolToSmiles(frag_mol)

        frag_smi_list = frag_smi.split(".")

        side_chain_smi_list = []

        for s in frag_smi_list:
            if s.count("*") == 1:
                replaced_string = re.sub(r"(\d+)\*", r"*:\1", s)  # [1*] -> [*:1]
                side_chain_smi_list.append(replaced_string)

        side_chain_smi = ".".join(side_chain_smi_list)

        side_chain = Chem.MolFromSmiles(side_chain_smi)

        yield from gen_new_replacements(
            [(env_smarts, core_smi, atom_ids)],
            side_chain,
            db_manager,
            min_inc,
            max_inc,
            max_replacements,
            num_cpus,
            radius=radius,
            dist=dist,
            min_freq=min_freq,
            products=products,
        )


def grow_mol(
    mol,
    db_manager,
    radius=3,
    min_inc=1,
    max_inc=2,
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
        ids = set(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1).difference(
            ids
        )  # ids of Hs to protect
        protected_ids.update(ids)  # Hs should be protected

    return mutate_mol(
        mol,
        db_manager,
        radius,
        min_inc=min_inc,
        max_inc=max_inc,
        max_replacements=max_replacements,
        replace_ids=None,
        protected_ids=protected_ids,
        num_cpus=num_cpus,
        symmetry_fixes=symmetry_fixes,
    )



def link_mol(
    mol1,
    mol2,
    db_manager,
    radius=3,
    min_inc=1,
    max_inc=2,
    dist=None,
    min_freq=0,
    max_replacements=None,
    replace_ids_1=None,
    replace_ids_2=None,
    symmetry_fixes=False,
    num_cpus=1,
    attach_id=None,
):
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
        ids = set(a.GetIdx() for a in mol1.GetAtoms() if a.GetAtomicNum() == 1).difference(
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
        ids = set(a.GetIdx() for a in mol2.GetAtoms() if a.GetAtomicNum() == 1).difference(
            ids
        )  # ids of Hs to protect
        protected_ids_2.update(ids)  # Hs should be protected

    fragments = fragment_mol_link(
        mol1, mol2, radius, protected_ids_1=protected_ids_1, protected_ids_2=protected_ids_2,
    )  # [(env smiles, core smiles, list of atom ids)]
    print("fragments", fragments)

    mol = Chem.CombineMols(mol1, mol2)
    mol_copy = deepcopy(mol)

    yield from gen_new_replacements(
        fragments,
        mol_copy,
        db_manager,
        min_inc,
        max_inc,
        max_replacements,
        num_cpus,
        radius,
        dist=dist,
        min_freq=min_freq,
    )

