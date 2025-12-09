#!/usr/bin/env python
# -*- encoding: utf-8 -*-


import copy
import logging
from multiprocessing import Pool, cpu_count

import numpy as np
from rdkit import Chem
from rdkit.Chem.Descriptors import (TPSA, MolLogP, MolWt, NumHAcceptors,
                                    NumHDonors)
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds

from ta_gen.tool.property_utils import GetNumHeavyAtoms
from ta_gen.utils.common_utils import fault_tolerance

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)


criteria_method = {
    "MW": MolWt,
    "LogP": MolLogP,
    "HBD": NumHDonors,
    "HBA": NumHAcceptors,
    "PSA": TPSA,
    "Rot": CalcNumRotatableBonds,
    "HAC": GetNumHeavyAtoms,
}


class Filter(object):

    @staticmethod
    def _compute_mol_properties(smi, properties):
        """
        Compute properties of input molecule

        Args:
            smi (str): smiles
            properties (list): properties

        Returns:
            list: properties of input molecule

        """
        mol = Chem.MolFromSmiles(smi)
        if mol:
            return tuple(criteria_method[col](mol) for col in properties)
        else:
            return (-999,) * len(properties)

    def compute_properties(self, df, properties, threads=1):
        """
        Compute properties

        Args:
            df (pandas.DataFrame): input data
            threads (int): number of threads

        Returns:
            df (pandas.DataFrame): output data

        """
        pool = Pool(threads)
        ret = []
        for _, smi in enumerate(df["smiles"]):
            ret.append(
                pool.apply_async(self._compute_mol_properties, (smi, properties))
            )
        pool.close()
        pool.join()

        df[properties] = np.asarray([x.get() for x in ret])
        return df

    def compute_violation_count(self, df, criteria_boundary):
        valid_cnt = 0
        for col, (low, high) in criteria_boundary.items():
            valid_cnt += df[col].between(low, high).astype(int)
        violation_cnt = len(criteria_boundary) - valid_cnt
        return violation_cnt

    def filter_by_properties(
        self, df, criteria_boundary, num_of_violations=0, keep_filter_columns=False
    ):
        """
        Filter by properties
        Args:
            df (pandas.DataFrame): input data
            criteria_boundary (dict): criteria
        Returns:
            df (pandas.DataFrame): output data
        """
        filter_columns = list(criteria_boundary.keys())
        df = self.compute_properties(df, filter_columns)
        violation_cnt = self.compute_violation_count(df, criteria_boundary)
        df = df[violation_cnt <= num_of_violations]
        if not keep_filter_columns:
            df.drop(filter_columns, axis=1, inplace=True)
        return df

    def filter_on_df(
        self,
        df,
        num_of_violations=0,
        threads=1,
        **kwargs,
    ):
        if df.shape[0] == 0:
            logging.info("Library size = 0 - Skip filtering")
            return df

        required_columns = ["MW", "LogP", "HBD", "HBA", "PSA", "Rot"]
        lacked_columns = [key for key in required_columns if key not in df.columns]
        print(f"lacked_columns: {lacked_columns or 'None'}")
        if lacked_columns:
            if threads <= 0:
                threads = cpu_count()
            logging.info(f"Computing properties from library on {threads} threads")
            df = self.compute_properties(df, lacked_columns, threads)
            logging.info("Computing properties from library")

        print(f"Library size before filtering: {len(df)}")
        logging.info("Filtering library based on properties")

        criteria_boundary = {
            "MW": (kwargs.get("mw_min", 0), kwargs.get("mw_max", 500)),
            "LogP": (kwargs.get("logP_min", -5), kwargs.get("logP_max", 5)),
            "HBD": (kwargs.get("hbd_min", 0), kwargs.get("hbd_max", 5)),
            "HBA": (kwargs.get("hba_min", 0), kwargs.get("hba_max", 10)),
            "PSA": (kwargs.get("tpsa_min", 0), kwargs.get("tpsa_max", 200)),
            "Rot": (kwargs.get("rot_min", 0), kwargs.get("rot_max", 10)),
        }

        filtered_df = self.filter_by_properties(
            df, criteria_boundary, num_of_violations
        )
        print(f"Library size after filtering: {filtered_df.shape[0]}")
        return filtered_df

    @staticmethod
    def add_to_argparse(parser):
        group_filter = parser.add_argument_group("filter options")
        group_filter.add_argument(
            "--hba-min",
            type=float,
            default=0,
            help="Min number of hydrogen bond acceptor",
        )
        group_filter.add_argument(
            "--hba-max",
            type=float,
            default=10,
            help="Max number of hydrogen bond acceptor",
        )
        group_filter.add_argument(
            "--hbd-min", type=float, default=0, help="Min number of hydrogen bond donor"
        )
        group_filter.add_argument(
            "--hbd-max", type=float, default=5, help="Max number of hydrogen bond donor"
        )
        group_filter.add_argument(
            "--mw-min", type=float, default=0, help="Min molecular weight"
        )
        group_filter.add_argument(
            "--mw-max", type=float, default=500, help="Max molecular weight"
        )
        group_filter.add_argument("--logP-min", type=float, default=-5, help="Min logP")
        group_filter.add_argument("--logP-max", type=float, default=5, help="Max logP")
        group_filter.add_argument("--tpsa-min", type=float, default=0, help="Min TPSA")
        group_filter.add_argument(
            "--tpsa-max", type=float, default=200, help="Max TPSA"
        )
        group_filter.add_argument(
            "--rot-min", type=float, default=0, help="Min rotatable bonds"
        )
        group_filter.add_argument(
            "--rot-max", type=float, default=10, help="Max rotatable bonds"
        )
        group_filter.add_argument(
            "--num_of_violations", type=int, default=0, help="Number of violations"
        )
        group_filter.add_argument(
            "--threads",
            type=int,
            default=1,
            help="Number of threads, 0 for all threads",
        )


class PreFilter(object):
    criteria_method = {}

    @staticmethod
    def _compute_mol_properties(smi, scaffold, properties):
        """
        Compute properties of input molecule

        Args:
            smi (str): smiles
            properties (list): properties

        Returns:
            list: properties of input molecule

        """
        mol = Chem.MolFromSmiles(smi)
        scaffold = Chem.MolFromSmiles(scaffold)
        if mol:
            return tuple(
                PreFilter.criteria_method[col](mol, scaffold) for col in properties
            )
        else:
            return (-999,) * len(properties)

    @fault_tolerance
    def get_MW_prefilter(self, reagent_smi, reagent_order, scaffold_mw):
        """
        Compute MW of input molecule

        Args:
            reagent_smi (str): reagent smiles
            reagent_order (dict): info of reagent
            scaffold_mw (float): MW of scaffold


        Returns:
            float: MW of input molecule

        """
        reagent = Chem.MolFromSmiles(reagent_smi)
        leaving_group_pattern = Chem.MolFromSmarts(
            reagent_order["reagent_leaving_group"]
        )
        reagent_mw = MolWt(Chem.DeleteSubstructs(reagent, leaving_group_pattern))

        return reagent_mw + scaffold_mw

    @fault_tolerance
    def get_HAC_prefilter(self, reagent_smi, reagent_order, scaffold_hac):
        """
        Compute HAC of input molecule

        Args:
            reagent_smi (str): reagent smiles
            reagent_order (dict): info of reagent
            scaffold_hac (int): HAC of scaffold

        Returns:
            int: HAC of input molecule

        """
        reagent = Chem.MolFromSmiles(reagent_smi)
        leaving_group_pattern = Chem.MolFromSmarts(
            reagent_order["reagent_leaving_group"]
        )
        reagent_hac = GetNumHeavyAtoms(
            Chem.DeleteSubstructs(reagent, leaving_group_pattern)
        )

        return reagent_hac + scaffold_hac

    def compute_properties(self, df, reagent_order, properties, threads=1):
        properties = copy.deepcopy(properties)
        with Pool(threads) as pool:
            if "MW" in properties:
                scaffold = Chem.MolFromSmiles(reagent_order["scaffold"])
                leaving_group_pattern = Chem.MolFromSmarts(
                    reagent_order["scaffold_leaving_group"]
                )
                scaffold_mw = MolWt(
                    Chem.DeleteSubstructs(scaffold, leaving_group_pattern)
                )
                ret = []
                for _, smi in enumerate(df["smiles"]):
                    ret.append(
                        pool.apply_async(
                            self.get_MW_prefilter, (smi, reagent_order, scaffold_mw)
                        )
                    )
                df["MW"] = np.asarray([x.get() for x in ret])
                properties.remove("MW")

            if "HAC" in properties:
                scaffold = Chem.MolFromSmiles(reagent_order["scaffold"])
                leaving_group_pattern = Chem.MolFromSmarts(
                    reagent_order["scaffold_leaving_group"]
                )
                scaffold_hac = GetNumHeavyAtoms(
                    Chem.DeleteSubstructs(scaffold, leaving_group_pattern)
                )

                ret = []
                for _, smi in enumerate(df["smiles"]):
                    ret.append(
                        pool.apply_async(
                            self.get_HAC_prefilter, (smi, reagent_order, scaffold_hac)
                        )
                    )
                df["HAC"] = np.asarray([x.get() for x in ret])
                properties.remove("HAC")

            # rest properties
            ret = []
            for _, smi in enumerate(df["smiles"]):
                ret.append(
                    pool.apply_async(
                        self._compute_mol_properties,
                        (smi, reagent_order["scaffold"], properties),
                    )
                )

            df[properties] = np.asarray([x.get() for x in ret])
        return df

    def compute_violation_count(self, df, criteria_boundary):
        valid_cnt = 0
        for col, (low, high) in criteria_boundary.items():
            valid_cnt += df[col].between(low, high).astype(int)
        violation_cnt = len(criteria_boundary) - valid_cnt
        return violation_cnt

    def filter_by_properties(
        self,
        df,
        reagent_order,
        criteria_boundary,
        num_of_violations=0,
        keep_filter_columns=False,
        sequential=True,
    ):
        """
        Filter by properties
        Args:
            df (pandas.DataFrame): input data
            reagent_order (dict): reagent order
            criteria_boundary (dict): criteria
        Returns:
            df (pandas.DataFrame): output data
        """
        filter_columns = list(criteria_boundary.keys())
        if sequential:
            for col, (low, high) in criteria_boundary.items():
                df = self.compute_properties(df, reagent_order, [col])
                df = df[df[col].between(low, high)]
        else:
            df = self.compute_properties(df, reagent_order, filter_columns)
            violation_cnt = self.compute_violation_count(df, criteria_boundary)
            df = df[violation_cnt <= num_of_violations]
        if not keep_filter_columns:
            df.drop(filter_columns, axis=1, inplace=True)
        return df
