#!/usr/bin/env python
# -*- coding:utf-8 -*-


from rdkit import Chem
from rdkit.Chem import rdMMPA
from ta_gen.crem_utils.mol_context import get_canon_context_core, patt_remove_map, patt_remove_brackets
from collections import defaultdict
from itertools import product


def get_atom_prop(molecule, prop="Index"):
    res = []
    for a in molecule.GetAtoms():
        if a.GetAtomicNum():
            res.append(a.GetIntProp(prop))
    return tuple(sorted(res))


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
            ids_list = [set(i) for i in mol.GetSubstructMatches(Chem.MolFromSmarts(smi))]
            for ids_matched in ids_list:
                for ids_eq in product(
                    *(atom_eq[i] for i in item[2])
                ):  # enumerate all combinations of equivalent atoms
                    if ids_matched == set(ids_eq):
                        extended_output.append((item[0], item[1], tuple(sorted(ids_eq))))
    return extended_output


def fragment_mol(  # noqa: C901
    mol, radius=3, return_ids=True, keep_stereo=False, protected_ids=None, symmetry_fixes=False
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


def __fragment_link_one_mol(mol, return_ids=True, protected_ids=None, keep_stereo=False):
    if return_ids:
        for atom in mol.GetAtoms():
            atom.SetIntProp("Index", atom.GetIdx())

    frags = rdMMPA.FragmentMol(
        mol, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100
    )

    if protected_ids:
        filtered_frags = []
        protected_ids = set(protected_ids)
        for _, chains in frags:
            for atom in chains.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    for d in atom.GetNeighbors():
                        if d.GetAtomicNum() != 1 and d.GetIdx() not in protected_ids:
                            filtered_frags.append((None, chains))
    else:
        filtered_frags = frags

    ls = []
    for _, chains in filtered_frags:
        ids = []
        for atom in chains.GetAtoms():
            if atom.GetAtomicNum() == 0:
                for d in atom.GetNeighbors():
                    if d.GetAtomicNum() == 1:
                        ids = [d.GetIntProp("Index")]
            if ids:
                break  # only one such occurrence can be
        a, b = Chem.MolToSmiles(chains, isomericSmiles=keep_stereo).split(".")
        if a == "[H][*:1]":
            ls.append([b, ids])
        else:
            ls.append([a, ids])
    return ls



def fragment_mol_link(
    mol1,
    mol2,
    radius=3,
    keep_stereo=False,
    protected_ids_1=None,
    protected_ids_2=None,
    return_ids=True,
):
    if protected_ids_1 or protected_ids_2:
        return_ids = True

    if return_ids:
        for atom in mol1.GetAtoms():
            atom.SetIntProp("Index", atom.GetIdx())
        for atom in mol2.GetAtoms():
            atom.SetIntProp("Index", atom.GetIdx())

    frags_1 = __fragment_link_one_mol(mol1, return_ids, protected_ids_1, keep_stereo)

    frags_2 = __fragment_link_one_mol(mol2, return_ids, protected_ids_2, keep_stereo)

    for i in range(len(frags_1)):
        frags_1[i][0] = frags_1[i][0].replace("*:1", "*:2")

    q = []
    for (fr1, ids1), (fr2, ids2) in product(frags_1, frags_2):
        q.append(["%s.%s" % (fr1, fr2), ids1, ids2])

    fake_core = "[*:1]C[*:2]"
    output = []

    for chains, ids_1, ids_2 in q:
        env, frag = get_canon_context_core(
            chains, fake_core, radius=radius, keep_stereo=keep_stereo
        )
        output.append((env, "[H][*:1].[H][*:2]", ids_1, ids_2))

    return output  # list of tuples (env smiles, core smiles, list of atom ids)