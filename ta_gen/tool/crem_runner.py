#!/usr/bin/env python
# -*- encoding: utf-8 -*-
import os
import re

from rdkit import Chem

from ta_gen.tool.crem import grow_mol2, update_db_config
from ta_gen.tool.rgroup_enumeration_tool import RGroupEnumerationTool
from ta_gen.utils.common_utils import standardize_smiles, validate_smiles
from ta_gen.utils.const import MAXINUM_NUM_OF_OUTPUT_MOLS
from ta_gen.utils.logger import LOGGER


class CremRunner(RGroupEnumerationTool):
    """
    This class runs the CReM for generating molecules.

    Inherits from `RGroupEnumerationTool`.

    Args:
        scaffold (str): scaffold string
        db_config (str): fragments db config path.
        num_crem_cycles (int, optional): Number of times to run the CReM algorithm.
            Defaults to 2.
        **kwargs: Keyword arguments to pass to `RGroupEnumerationTool`.

    """

    def __init__(self, scaffold: str, db_config: str, db_engine: str, **kwargs):
        super().__init__(scaffold=scaffold)
        self.db_config = update_db_config(db_engine, db_config)
        self.radius = kwargs.get("radius", 3)
        self.min_atoms = kwargs.get("min_atoms", 1)
        self.max_atoms = kwargs.get("max_atoms", 2)
        self.result_out = kwargs.get("result_out", "result.csv")
        self.num_cpus = kwargs.get("num_cpus", 1)
        self.work_dir = kwargs.get("work_dir", os.getcwd())
        self.max_mols_gen = kwargs.get("max_mols_gen", MAXINUM_NUM_OF_OUTPUT_MOLS)
        self.db_query = kwargs.get("db_query", {})

    @staticmethod
    def validate_product(initial_smart, smiles):
        if not validate_smiles(smiles):
            return False
        mol = Chem.MolFromSmiles(smiles)
        if not mol.HasSubstructMatch(initial_smart):
            return False
        return True

    def additional_grow(self, smi, smarts, protect_id, max_replacements):
        out_smis = set()
        pre_out_mol = Chem.MolFromSmiles(smi)
        pre_out_mol = Chem.AddHs(pre_out_mol)
        subs = pre_out_mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
        for j in range(len(subs)):
            replace_id = subs[j][protect_id]
            print("replace id", replace_id)
            smis = grow_mol2(
                pre_out_mol,
                self.db_config,
                replace_ids=[replace_id],
                min_atoms=self.min_atoms,
                max_atoms=self.max_atoms,
                radius=self.radius,
                symmetry_fixes=True,
                max_replacements=max_replacements,
                **self.db_query,
            )
            out_smis |= set(smis)
        return out_smis

    def format_attachment_point_nums_in_scaffold_str(self, scaffold):
        pattern = r"\[\*\:(\d+)\]"
        replacement = "[*:{new_idx}]"
        scaffold = re.sub(
            pattern,
            lambda match: replacement.format(new_idx=int(match.group(1)) + 1),
            scaffold,
        )
        print(f"scaffold attach index starts from 1: {scaffold}")

        matches = re.findall(pattern, scaffold)
        matches = [int(i) for i in matches]
        print("attachment points", matches)  # [1, 2, 3]
        return scaffold, matches

    def convert_hydrogen_to_star(self, mol):
        mol_h = Chem.AddHs(mol)
        for atom in mol_h.GetAtoms():
            if atom.GetAtomicNum() == 1 and atom.GetAtomMapNum() == 0:
                atom.SetAtomicNum(0)
        smi = Chem.MolToSmiles(mol_h)
        print("convert H to *", smi)
        mol_h = Chem.MolFromSmiles(smi)
        return mol_h

    def get_attach_id_list(self, mol_h, matches):
        attach_id_list = []
        for atom in mol_h.GetAtoms():
            if atom.GetAtomicNum() != 0:  # skip H
                continue
            if atom.GetAtomMapNum() in matches:
                attach_id_list.append((atom.GetAtomMapNum(), atom.GetIdx()))

        attach_id_list.sort(key=lambda x: x[0])
        attach_id_list = [x[1] for x in attach_id_list]
        print("attachment atom index", attach_id_list)
        return attach_id_list

    def get_protect_id_list(self, attach_id_list, mol_h):
        protect_id_list = []
        for i in attach_id_list:
            neighbors = [x.GetIdx() for x in mol_h.GetAtomWithIdx(i).GetNeighbors()]
            protect_id_list.extend(neighbors)
        print("protect atom index", protect_id_list)
        return protect_id_list

    def convert_star_to_hydrogen(self, mol_h):
        for atom in mol_h.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomicNum(1)
                atom.SetAtomMapNum(0)

    def prepare(self):
        scaffold, matches = self.format_attachment_point_nums_in_scaffold_str(
            self.scaffold
        )

        mol = Chem.MolFromSmiles(standardize_smiles(scaffold))

        # convert H to *
        mol_h = self.convert_hydrogen_to_star(mol)

        # get attchement point id
        attach_id_list = self.get_attach_id_list(mol_h, matches)

        # get neighbors
        protect_id_list = self.get_protect_id_list(attach_id_list, mol_h)

        # convert * to H
        self.convert_star_to_hydrogen(mol_h)

        return mol_h, protect_id_list, attach_id_list

    def grow(self, mol_h, protect_id_list, attach_id_list, max_replacements):
        # init grow
        LOGGER.info("Start initial grow...")
        LOGGER.info(f"replace_id: {protect_id_list[0]}")
        LOGGER.info(f"max_replacements: {max_replacements}")
        LOGGER.info(f"db_query: {self.db_query}")
        out_smis = grow_mol2(
            mol_h,
            self.db_config,
            replace_ids=[protect_id_list[0]],
            min_atoms=self.min_atoms,
            max_atoms=self.max_atoms,
            radius=self.radius,
            symmetry_fixes=True,
            ncores=self.num_cpus,
            max_replacements=max_replacements,
            attach_id=attach_id_list[0],
            **self.db_query,
        )
        LOGGER.info("Initial grow done")
        LOGGER.info(f"attatchment point {attach_id_list[0]} output: {len(out_smis)}")
        return out_smis

    @staticmethod
    def add_to_argparse(parser):
        """
        The add_to_argparse function adds arguments related to the Crem module to an
        argparse.ArgumentParser object.

        Args:
            parser (argparse.ArgumentParser): An instance of argparse.ArgumentParser to add
            arguments to.

        Returns:
            argparse.ArgumentParser: The same argparse.ArgumentParser object with additional
            arguments added.
        """
        group_crem = parser.add_argument_group("crem options")
        group_crem.add_argument("--db_config", type=str, required=True)
        group_crem.add_argument("--db_engine", type=str, default="postgres")
        group_crem.add_argument("--radius", type=int, default=3)
        group_crem.add_argument("--min_atoms", type=int, default=1)
        group_crem.add_argument("--max_atoms", type=int, default=2)
        group_crem.add_argument("--num_cpus", type=int, default=1)
        group_crem.add_argument("--db_query", type=str)
        return parser

    def clean_tmp_files(self):
        # No files to clean
        pass
