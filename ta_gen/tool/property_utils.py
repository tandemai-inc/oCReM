#!/usr/bin/env python
# -*- coding:utf-8 -*-


from rdkit.Chem.Descriptors import MolWt


def GetNumHeavyAtoms(mol):
    return mol.GetNumHeavyAtoms()


def get_MW_prefilter(scaffold, reactant, reaction):
    """
    Compute MW of scaffold + reactant + reaction
    Args:
        scaffold (rdkit.Chem.rdchem.Mol): scaffold
        reactant (rdkit.Chem.rdchem.Mol): reactant
        reaction (rdkit.Chem.rdchem.Mol): reaction
    Returns:
        float: MW of scaffold + reactant
    """

    return MolWt(scaffold) + MolWt(reactant) + MolWt(reaction)


def get_HAC_prefilter(scaffold, reactant):
    return GetNumHeavyAtoms(scaffold) + GetNumHeavyAtoms(reactant)
