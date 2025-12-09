#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import itertools
import math
from functools import partial
from multiprocessing import Pool

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

from ta_gen.utils.common_utils import (
    calculate_mol_properties,
    get_random_sample_from_df,
    load_ionizable_structure_smarts,
)
from ta_gen.utils.const import DISPLAY_ORDER
from ta_gen.utils.logger import LOGGER
from ta_gen.utils.prop_filter import Filter
from ta_gen.utils.rgroup_identification import identify_rgroups
from ta_gen.utils.structural_alerts import (
    generate_tafilter_id_by_count,
    get_scaffold_pattern_count,
    get_sdf_tafilter_id,
    get_side_chain_pattern_count,
)

RDLogger.DisableLog("rdApp.*")
CHUNK_SIZE = 500000


def chunkify(lst, n):
    """Split list into n-sized chunks."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def _get_properties_and_tafilter_id(
    smi,
    cnt_scaffold_pattern,
    cnt_side_chain_pattern,
    ignore_tf_alerts,
    ionizable_structure_smarts,
):
    mol_properties = {"smiles": smi}
    mol_properties.update(calculate_mol_properties(smi, ionizable_structure_smarts))
    tafilter_id = generate_tafilter_id_by_count(
        smi, cnt_scaffold_pattern, cnt_side_chain_pattern, ignore_tf_alerts
    )
    mol_properties["tafilter_id"] = ";".join(tafilter_id)
    return mol_properties


class PostProcessor(object):

    @staticmethod
    def get_properties_and_tafilter_id(
        result_df,
        scaffold,
        num_cpus,
        side_chain_smi,
        ref_hac,
        ignore_tf_alerts,
        chunk_size=CHUNK_SIZE,
    ):
        if result_df.empty:
            return result_df
        # Calculate properties
        LOGGER.info("Start calculating properties")
        smi_list = result_df["smiles"].to_list()

        # Filter scaffold or side chain
        cnt_scaffold_pattern = get_scaffold_pattern_count(scaffold)
        cnt_side_chain_pattern = get_side_chain_pattern_count(side_chain_smi)

        LOGGER.info(f"chunk size: {chunk_size}")
        LOGGER.info("Start calculating properties on chunks")
        ionizable_structure_smarts = load_ionizable_structure_smarts()
        __get_properties_and_tafilter_id = partial(
            _get_properties_and_tafilter_id,
            cnt_scaffold_pattern=cnt_scaffold_pattern,
            cnt_side_chain_pattern=cnt_side_chain_pattern,
            ignore_tf_alerts=ignore_tf_alerts,
            ionizable_structure_smarts=ionizable_structure_smarts,
        )
        with Pool(num_cpus) as pool:
            res_list = []
            # Avoid memory leak
            for chunk in chunkify(smi_list, chunk_size):
                res_list.extend(pool.map(__get_properties_and_tafilter_id, chunk))

        properties_df = pd.DataFrame(res_list)
        result_df = result_df.merge(properties_df, on="smiles", suffixes=("", "_new"))
        # Calculate HAC diff
        if ref_hac:
            result_df["HAC_diff"] = result_df["numA"] - ref_hac
        else:
            result_df["HAC_diff"] = 0

        return result_df

    @staticmethod
    def get_rgroups(
        smi_list, scaffold, attachment_points_mapping, num_cpus, chunk_size=CHUNK_SIZE
    ):
        LOGGER.info("Start identifying R-groups")
        num_chunks = math.ceil(len(smi_list) / chunk_size)
        LOGGER.info(f"total chunks: {num_chunks}")
        updated_smi_list = []
        rgroups_list = []
        for i in range(num_chunks):
            LOGGER.info(f"Processing chunk {i + 1}")
            chunk = smi_list[i * chunk_size : (i + 1) * chunk_size]
            chunk_smis, chunk_rgroups = identify_rgroups(
                chunk, [scaffold], attachment_points_mapping, num_cpus=num_cpus
            )
            updated_smi_list.extend(chunk_smis)
            rgroups_list.extend(chunk_rgroups)

        return updated_smi_list, rgroups_list

    @staticmethod
    def post_process(
        result_data,
        run_ta_prop,
        tool_name,
        attachment_points_mapping,
        num_cpus,
        max_mols_gen,
        result_file,
        ta_prop_args=None,
        scaffold=None,
        side_chain_smi=None,
        ref_sdf=None,
        remove_radical=True,
    ):
        if max_mols_gen <= 0:
            max_mols_gen = math.inf

        column_order = DISPLAY_ORDER[tool_name]

        result_df = PostProcessor.initialize_result_df(
            result_data, column_order, result_file
        )
        if result_df is None:
            return

        result_df.drop_duplicates(subset="smiles", keep="first", inplace=True)

        if run_ta_prop:
            result_df = PostProcessor.apply_filter(
                result_df, ta_prop_args, column_order, result_file
            )
            if result_df is None:
                return

        if remove_radical:
            result_df = PostProcessor.remove_radicals(result_df)

        result_df = get_random_sample_from_df(result_df, max_mols_gen)

        ignore_tf_alerts, ref_hac = PostProcessor.get_ignore_tf_alerts_and_ref_hac(
            ref_sdf, scaffold
        )

        if tool_name.upper() in ["CREM", "LIB-INVENT"]:
            result_df = PostProcessor.identify_rgroups(
                result_df, scaffold, attachment_points_mapping, num_cpus
            )
        elif tool_name.upper() == "CORE-HOPPING":
            result_df = PostProcessor.get_ch_core_properties(result_df, num_cpus)

        result_df = PostProcessor.get_properties_and_tafilter_id(
            result_df, scaffold, num_cpus, side_chain_smi, ref_hac, ignore_tf_alerts
        )

        if result_df.empty:
            raise ValueError("No molecules generated.")

        PostProcessor.write_output(result_df, result_file, column_order, tool_name)

    @staticmethod
    def initialize_result_df(result_data, column_order, result_file):
        if isinstance(result_data, pd.DataFrame):
            result_df = result_data
        else:
            result_df = pd.DataFrame(result_data)

        if len(result_df) == 0:
            LOGGER.warning("No results to process since no molecules generated.")
            with open(result_file, "w") as f:
                f.write(",".join(column_order) + "\n")
            return None
        return result_df

    @staticmethod
    def apply_filter(result_df, ta_prop_args, column_order, result_file):
        result_df = Filter().filter_on_df(result_df, **ta_prop_args)
        if len(result_df) == 0:
            LOGGER.warning("All molecules are filtered out.")
            with open(result_file, "w") as f:
                f.write(",".join(column_order) + "\n")
            return None
        return result_df

    @staticmethod
    def remove_radicals(result_df):
        LOGGER.info("Remove radical molecules")
        mask = result_df["smiles"].apply(
            lambda x: Descriptors.NumRadicalElectrons(Chem.MolFromSmiles(x)) == 0
        )
        return result_df[mask]

    @staticmethod
    def get_ignore_tf_alerts_and_ref_hac(ref_sdf, scaffold):
        ignore_tf_alerts = set()
        if ref_sdf:
            ref_hac = Chem.MolFromMolFile(ref_sdf).GetNumHeavyAtoms()
            ignore_tf_alerts = get_sdf_tafilter_id(ref_sdf)
            LOGGER.info("Ignore TF alerts: " + ",".join(ignore_tf_alerts))
        elif scaffold:
            ref_hac = Chem.MolFromSmiles(scaffold).GetNumHeavyAtoms()
        else:
            ref_hac = 0  # Skip
        return ignore_tf_alerts, ref_hac

    @staticmethod
    def identify_rgroups(result_df, scaffold, attachment_points_mapping, num_cpus):
        smi_list = result_df["smiles"].to_list()
        updated_smi_list, rgroups_list = PostProcessor.get_rgroups(
            smi_list, scaffold, attachment_points_mapping, num_cpus
        )
        return pd.DataFrame({"smiles": updated_smi_list, "rgroups": rgroups_list})

    @staticmethod
    def write_output(result_df, result_file, column_order, tool_name):
        LOGGER.info(f"Write output to {result_file}")

        if "mol_id" not in result_df.columns:
            result_df["mol_id"] = range(1, len(result_df) + 1)

        if tool_name != "post-process":
            result_df = result_df.reindex(columns=column_order)

        result_df.to_csv(result_file, index=False)
        LOGGER.info(f"Complete! Total output: {result_df.shape[0]}")

    @staticmethod
    def effective_length(mol: Chem.Mol) -> int:
        # if a single atom has more than one attachment point, return 0
        attachement_indices = []
        num_attachments = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                num_attachments += 1
                attachement_indices.append(atom.GetIdx())
        if len(attachement_indices) < num_attachments:
            return 0
        else:
            pairs = itertools.combinations(attachement_indices, 2)
            return int(min(len(Chem.GetShortestPath(mol, i, j)) - 1 for i, j in pairs))

    @staticmethod
    def get_ch_core_properties(result_df, num_cpus, chunk_size=CHUNK_SIZE):
        core_smi_list = result_df["core"].unique().tolist()
        LOGGER.info("Start calculating core properties")
        LOGGER.info(f"chunk size: {chunk_size}")
        LOGGER.info("Start calculating properties on chunks")
        with Pool(num_cpus) as pool:
            properties_list = list(
                pool.map(
                    PostProcessor.get_new_core_data, core_smi_list, chunksize=chunk_size
                )
            )
            properties_df = pd.DataFrame(properties_list)
        result_df = result_df.merge(
            properties_df, on="core", suffixes=("", "_new"), how="left"
        )
        return result_df

    @staticmethod
    def get_new_core_data(core_smi):
        # core_smi: [*:2]CC1(COCCC)CCCC1[*:1]
        # get Core Effective Length, Core Num. Aro. Rings, Core Num. Heavy Atoms
        core = Chem.MolFromSmiles(core_smi)
        effective_length = PostProcessor.effective_length(core)
        num_rings = Descriptors.NumAromaticRings(core)
        num_ha = core.GetNumHeavyAtoms()
        return {
            "core": core_smi,
            "coreEffectiveLength": effective_length,
            "coreNumAromaticRings": num_rings,
            "coreNumHeavyAtoms": num_ha,
        }
