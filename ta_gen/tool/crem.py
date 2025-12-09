# Pavel Polishchuk, 2017

import random
import re
import sqlite3
from collections import defaultdict
from copy import deepcopy
from itertools import product
from multiprocessing import Pool, cpu_count

import psycopg2

from ta_gen.utils.common_utils import get_db_config, standardize_smiles

from crem.mol_context import get_canon_context_core, patt_remove_map

from rdkit import Chem
from rdkit.Chem import rdMMPA

cycle_pattern = re.compile(r"[a-zA-Z\]][1-9]+")
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
patt_remove_brackets = re.compile(r"\(\)")


def __extend_output_by_equivalent_atoms(mol, output):
    atom_ranks = list(
        Chem.CanonicalRankAtoms(
            mol, breakTies=False, includeChirality=False, includeIsotopes=False
        )
    )
    tmp = defaultdict(list)
    for i, rank in enumerate(atom_ranks):
        tmp[rank].append(i)
    atom_eq = dict()  # dict of equivalent atoms
    for ids in tmp.values():
        if len(ids) > 1:
            for i in ids:
                atom_eq[i] = [j for j in ids if j != i]

    extended_output = []
    for item in output:
        if all(
            i in atom_eq.keys() for i in item[2]
        ):  # if all atoms of a fragment have equivalent atoms
            smi = patt_remove_map.sub("", item[1])
            smi = patt_remove_brackets.sub("", smi)
            ids_list = [
                set(i) for i in mol.GetSubstructMatches(Chem.MolFromSmarts(smi))
            ]
            for ids_matched in ids_list:
                for ids_eq in product(
                    *(atom_eq[i] for i in item[2])
                ):  # enumerate all combinations of equivalent atoms
                    if ids_matched == set(ids_eq):
                        extended_output.append(
                            (item[0], item[1], tuple(sorted(ids_eq)))
                        )
    return extended_output


def get_atom_prop(molecule, prop="Index"):
    res = []
    for a in molecule.GetAtoms():
        if a.GetAtomicNum():
            res.append(a.GetIntProp(prop))
    return tuple(sorted(res))


def __fragment_mol(  # noqa: C901
    mol,
    radius=3,
    return_ids=True,
    keep_stereo=False,
    protected_ids=None,
    symmetry_fixes=False,
):
    if protected_ids:
        return_ids = True

    # due to the bug https://github.com/rdkit/rdkit/issues/3040
    # outputs of rdMMPA.FragmentMol calls will contain duplicated fragments
    # they are removed by using this set
    output = set()

    # set original atom idx to keep them in fragmented mol
    if return_ids:
        for atom in mol.GetAtoms():
            atom.SetIntProp("Index", atom.GetIdx())

    # heavy atoms
    frags = rdMMPA.FragmentMol(
        mol, pattern="[!#1]!@!=!#[!#1]", maxCuts=4, resultsAsMols=True, maxCutBonds=30
    )
    frags += rdMMPA.FragmentMol(
        mol, pattern="[!#1]!@!=!#[!#1]", maxCuts=3, resultsAsMols=True, maxCutBonds=30
    )
    # hydrogen atoms
    frags += rdMMPA.FragmentMol(
        mol, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100
    )

    for i, (core, chains) in enumerate(frags):
        if core is None:  # single cut
            components = list(Chem.GetMolFrags(chains, asMols=True))
            ids_0 = get_atom_prop(components[0]) if return_ids else tuple()
            ids_1 = get_atom_prop(components[1]) if return_ids else tuple()
            if Chem.MolToSmiles(components[0]) != "[H][*:1]":  # context cannot be H
                env, frag = get_canon_context_core(
                    components[0], components[1], radius, keep_stereo
                )
                output.add((env, frag, ids_1))
            if Chem.MolToSmiles(components[1]) != "[H][*:1]":  # context cannot be H
                env, frag = get_canon_context_core(
                    components[1], components[0], radius, keep_stereo
                )
                output.add((env, frag, ids_0))
        else:  # multiple cuts
            # there are no checks for H needed because H can be present only in single cuts
            env, frag = get_canon_context_core(chains, core, radius, keep_stereo)
            output.add((env, frag, get_atom_prop(core) if return_ids else tuple()))

    if symmetry_fixes:
        extended_output = __extend_output_by_equivalent_atoms(mol, output)
        if extended_output:
            output.update(extended_output)

    if protected_ids:
        protected_ids = set(protected_ids)
        output = [item for item in output if protected_ids.isdisjoint(item[2])]

    return list(output)  # list of tuples (env smiles, core smiles, list of atom ids)


def __get_replacements_rowids_sqlite(db_cur, env, dist, min_atoms, max_atoms, radius):
    sql = (
        f"SELECT rowid FROM radius{radius} WHERE env = '{env}' "
        f"AND core_num_atoms BETWEEN {min_atoms} AND {max_atoms}"
    )
    if isinstance(dist, int):
        sql += f" AND fragment.dist2 = {dist}"
    elif isinstance(dist, tuple) and len(dist) == 2:
        sql += f" AND fragment.dist2 BETWEEN {dist[0]} AND {dist[1]}"

    db_cur.execute(sql)

    recs = db_cur.fetchall()

    return set(i[0] for i in recs)


def __get_replacements_rowids_pg(db_cur, env, dist, min_atoms, max_atoms):
    sql = (
        f"SELECT fragment.id "
        f"FROM env "
        "JOIN env_fragment ON env.id = env_fragment.env_id "
        "JOIN fragment ON fragment.id = env_fragment.fragment_id "
        f"WHERE env.name = '{env}' "
        f"AND fragment.core_num_atoms BETWEEN {min_atoms} AND {max_atoms} "
    )
    if isinstance(dist, int):
        sql += f" AND fragment.dist2 = {dist}"
    elif isinstance(dist, tuple) and len(dist) == 2:
        sql += f" AND fragment.dist2 BETWEEN {dist[0]} AND {dist[1]}"

    db_cur.execute(sql)

    recs = db_cur.fetchall()

    return set(i[0] for i in recs)


def __get_replacements(db_cur, row_ids, radius=2, engine="postgres", batch_size=100000):
    results = []
    for i in range(0, len(row_ids), batch_size):
        batch = row_ids[i : i + batch_size]
        if engine == "postgres":
            placeholders = ",".join("%s" for _ in batch)
            sql = f"SELECT id, core_smi FROM fragment WHERE id IN ({placeholders})"
        else:
            placeholders = ",".join("?" for _ in batch)
            sql = f"SELECT rowid, core_smi FROM radius{radius} WHERE rowid IN ({placeholders})"

        db_cur.execute(sql, batch)
        results.extend(db_cur.fetchall())
    return results


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


def gen_new_replacements(  # noqa: C901
    f,
    mol,
    db_config,
    min_inc,
    max_inc,
    max_replacements,
    ncores,
    radius,
    dist=None,
    products=None,
    return_core=False,
):
    if not products:
        products = set()

    if ncores == 1:
        for env_smarts, core_smi, _ in f:
            for new_core_smi in get_core_smi_replacements(
                db_config,
                env_smarts,
                core_smi,
                dist,
                min_inc,
                max_inc,
                max_replacements,
                radius,
            ):
                smi, _ = zip_new_replacement(new_core_smi, mol)
                if smi and smi not in products:
                    products.add(smi)
                    if return_core:
                        yield smi, new_core_smi
                    else:
                        yield smi
    else:
        with Pool(min(ncores, cpu_count())) as p:
            for env_smarts, core_smi, _ in f:
                for smi, new_core_smi in p.starmap(
                    zip_new_replacement,
                    (
                        (new_core_smi, mol)
                        for new_core_smi in get_core_smi_replacements(
                            db_config,
                            env_smarts,
                            core_smi,
                            dist,
                            min_inc,
                            max_inc,
                            max_replacements,
                            radius,
                        )
                    ),
                    chunksize=100,
                ):
                    if smi and smi not in products:
                        products.add(smi)
                        if return_core:
                            yield smi, new_core_smi
                        else:
                            yield smi


def mutate_mol(  # noqa: C901
    mol,
    db_config,
    radius=3,
    min_inc=-2,
    max_inc=2,
    max_replacements=None,
    replace_ids=None,
    protected_ids=None,
    symmetry_fixes=False,
    ncores=1,
    attach_id=None,
    attach_neighbor_id=None,
):
    protected_ids = update_protected_ids(mol, protected_ids, replace_ids)
    print("protected_ids", protected_ids)
    print("radius", radius)
    f = __fragment_mol(
        mol, radius, protected_ids=protected_ids, symmetry_fixes=symmetry_fixes
    )  # [(env smiles, core smiles, list of atom ids)]
    print("f", f)

    valid_f = []
    for t in f:
        num_heavy_atoms = Chem.MolFromSmiles(t[1]).GetNumHeavyAtoms()
        if num_heavy_atoms == 0:
            valid_f.append(t)

    print("valid f", f)

    mol_copy = deepcopy(mol)
    if attach_id:
        atom = mol_copy.GetAtomWithIdx(attach_id)
        if atom.GetAtomicNum() == 1:
            atom.SetAtomicNum(0)
            atom.SetAtomMapNum(1)
    elif attach_neighbor_id:
        for a in mol_copy.GetAtomWithIdx(attach_neighbor_id).GetNeighbors():
            if a.GetAtomicNum() == 1:  # set H atom to *1
                a.SetAtomicNum(0)
                a.SetAtomMapNum(1)
                break

    yield from gen_new_replacements(
        valid_f,
        mol_copy,
        db_config,
        min_inc,
        max_inc,
        max_replacements,
        ncores,
        radius,
    )


def get_core_smi_replacements(
    db_config, env_smarts, core_smi, dist, min_inc, max_inc, max_replacements, radius
):
    num_heavy_atoms = Chem.MolFromSmiles(core_smi).GetNumHeavyAtoms()
    min_atoms = num_heavy_atoms + min_inc
    max_atoms = num_heavy_atoms + max_inc

    engine = db_config["engine"]

    if engine == "postgres":
        pg_conn = psycopg2.connect(
            dbname=db_config.get("dbname"),
            user=db_config.get("user"),
            password=db_config.get("password"),
            host=db_config.get("host"),
            port=db_config.get("port"),
        )
        cur = pg_conn.cursor()

        row_ids = __get_replacements_rowids_pg(
            cur,
            env_smarts,
            dist,
            min_atoms,
            max_atoms,
        )
    elif engine == "sqlite":
        db_path = db_config["path"]
        sql_conn = sqlite3.connect(db_path)
        cur = sql_conn.cursor()
        row_ids = __get_replacements_rowids_sqlite(
            cur,
            env_smarts,
            dist,
            min_atoms,
            max_atoms,
            radius,
        )
    else:
        raise ValueError(f"Unsupported database engine: {engine}")

    if max_replacements is None or len(row_ids) <= max_replacements:
        row_ids = list(row_ids)
        res = __get_replacements(cur, row_ids, radius, engine)
    else:
        selected_row_ids = random.sample(row_ids, max_replacements)
        res = __get_replacements(cur, selected_row_ids, radius, engine)

    for _, new_core_smi in res:
        if new_core_smi != core_smi:
            yield new_core_smi


def zip_new_replacement(new_replacement, input_structure):
    try:
        new_replacement_mol = Chem.MolFromSmiles(new_replacement)
        final_mol = Chem.molzip(input_structure, new_replacement_mol)
        Chem.SanitizeMol(
            final_mol, catchErrors=True
        )  # Bonds can be restored to the aromatic bond type
        smi = Chem.MolToSmiles(final_mol, isomericSmiles=True)
        smi = standardize_smiles(smi)
        return smi, new_replacement
    except Exception:
        return None, None


def mol_to_smarts(mol, keep_h=True):
    # e.g. [H]-[CH2]-[*] -> [H]-[CH3]-[*]

    mol = Chem.Mol(mol)
    mol.UpdatePropertyCache()

    # change the isotope to 42
    for atom in mol.GetAtoms():
        if keep_h:
            s = sum(na.GetAtomicNum() == 1 for na in atom.GetNeighbors())
            if s:
                atom.SetNumExplicitHs(atom.GetTotalNumHs() + s)
        atom.SetIsotope(42)

    # print out the smiles - all the atom attributes will be fully specified
    smarts = Chem.MolToSmiles(mol, isomericSmiles=True, allBondsExplicit=True)
    # remove the 42 isotope labels
    smarts = re.sub(r"\[42", "[", smarts)

    return smarts


def combine_core_env_to_rxn_smarts(core, env, keep_h=True):  # noqa: C901
    if isinstance(env, str):
        m_env = Chem.MolFromSmiles(env, sanitize=False)
    if isinstance(core, str):
        m_frag = Chem.MolFromSmiles(core, sanitize=False)

    backup_atom_map = "backupAtomMap"

    # put all atom maps to atom property and remove them
    for a in m_env.GetAtoms():
        atom_map = a.GetAtomMapNum()
        if atom_map:
            a.SetIntProp(backup_atom_map, atom_map)
            a.SetAtomMapNum(0)
    for a in m_frag.GetAtoms():
        atom_map = a.GetAtomMapNum()
        if atom_map:
            a.SetIntProp(backup_atom_map, atom_map)
            a.SetAtomMapNum(0)

    # set canonical ranks for atoms in env without maps
    m_env.UpdatePropertyCache()
    for atom_id, rank in zip(
        [a.GetIdx() for a in m_env.GetAtoms()], list(Chem.CanonicalRankAtoms(m_env))
    ):
        a = m_env.GetAtomWithIdx(atom_id)
        if not a.HasProp(backup_atom_map):
            a.SetAtomMapNum(rank + 1)  # because ranks start from 0

    m = Chem.RWMol(Chem.CombineMols(m_frag, m_env))

    links = defaultdict(list)  # pairs of atom ids to create bonds
    att_to_remove = []  # ids of att points to remove
    for a in m.GetAtoms():
        if a.HasProp(backup_atom_map):
            i = a.GetIntProp(backup_atom_map)
            links[i].append(a.GetNeighbors()[0].GetIdx())
            att_to_remove.append(a.GetIdx())

    for i, j in links.values():
        m.AddBond(i, j, Chem.BondType.SINGLE)

    for i in sorted(att_to_remove, reverse=True):
        m.RemoveAtom(i)

    comb_sma = mol_to_smarts(m, keep_h)

    patt_remove_h = re.compile(r"(?<!\[H)H[1-9]*(?=:[0-9])")

    if not keep_h:  # remove H only in mapped env part
        comb_sma = patt_remove_h.sub("", comb_sma)
    return comb_sma


def mutate_mol_ch(
    mol,
    db_config,
    radius=3,
    min_inc=-2,
    max_inc=2,
    max_replacements=None,
    replace_ids=None,
    protected_ids=None,
    symmetry_fixes=False,
    ncores=1,
    side_chain=None,
):
    products = {Chem.MolToSmiles(mol)}

    protected_ids = update_protected_ids(mol, protected_ids, replace_ids)

    f = __fragment_mol(
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
                        mol.GetBondBetweenAtoms(
                            core_match[v[0]], core_match[v[1]]
                        ).GetIdx()
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
            db_config,
            min_inc,
            max_inc,
            max_replacements,
            ncores,
            radius=radius,
            products=products,
            return_core=True,
        )


def grow_mol(
    mol,
    db_config,
    radius=3,
    min_atoms=1,
    max_atoms=2,
    max_replacements=None,
    replace_ids=None,
    symmetry_fixes=False,
    ncores=1,
    attach_id=None,
):
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

    return mutate_mol(
        mol,
        db_config,
        radius,
        min_inc=min_atoms,
        max_inc=max_atoms,
        max_replacements=max_replacements,
        replace_ids=None,
        protected_ids=protected_ids,
        ncores=ncores,
        symmetry_fixes=symmetry_fixes,
        attach_id=attach_id,
        attach_neighbor_id=replace_ids[0] if replace_ids else None,
    )


def mutate_mol2(*args, **kwargs):
    """
    Convenience function which can be used to process molecules in parallel
    using multiprocessing module.
    It calls mutate_mol which cannot be used directly in multiprocessing because it is a generator

    :param args: positional arguments, the same as in mutate_mol function
    :param kwargs: keyword arguments, the same as in mutate_mol function
    :return: list with output molecules

    """
    return list(mutate_mol(*args, **kwargs))


def grow_mol2(*args, **kwargs):
    """
    Convenience function which can be used to process molecules in parallel
    using multiprocessing module.
    It calls grow_mol which cannot be used directly in multiprocessing because it is a generator

    :param args: positional arguments, the same as in grow_mol function
    :param kwargs: keyword arguments, the same as in grow_mol function
    :return: list with output molecules

    """
    return list(grow_mol(*args, **kwargs))


def update_db_config(db_engine, db_config):
    if db_engine == "postgres":
        db_config = get_db_config(db_config)
        db_config["engine"] = db_engine
    elif db_engine == "sqlite":
        db_config = {"path": db_config, "engine": db_engine}
    else:
        raise ValueError(f"Unsupported db_engine: {db_engine}")
    return db_config
