#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import re
from abc import ABC, abstractmethod

from ta_gen.utils.common_utils import validate_scaffold


class RGroupEnumerationTool(ABC):
    def __init__(self, scaffold: str):
        """
        Initializes an EnumerationTool object with scaffold string.

        Args:
            scaffold (str): scaffold string
        """
        self.attachment_points_mapping = {}
        self.scaffold = self.set_scaffold(scaffold)

    def set_scaffold(self, scaffold):
        """
        Setter for the scaffold attribute.

        Validates each scaffold in the input list and raises an exception if any are invalid.

        Args:
            scaffold (str): scaffold string
        """
        print("input scaffold", scaffold)
        if not validate_scaffold(scaffold):
            raise Exception(f"scaffold {scaffold} is not valid")
        formatted_scaffold, new2old, old2new = self.format_scaffold(scaffold)
        self.attachment_points_mapping[formatted_scaffold] = {
            "new2old": new2old,
            "old2new": old2new,
        }
        scaffold = formatted_scaffold
        print("updated scaffold", scaffold)
        print("attachment_points_mapping", self.attachment_points_mapping)
        return scaffold

    @staticmethod
    def format_scaffold(scaffold):
        """Make annotation indexes start from 0 and record the mapping relationship
        Args:
            scaffold (str): annotated smiles
        Returns:
            tuple: updated scaffold, mapping of new indexes and old indexes
        """
        #    [*:3]N1C(=O)NC2(C1=O)CC(C2)[*:6]
        # -> [*:0]N1C(=O)NC2(C1=O)CC(C2)[*:1]
        pattern = r"\[\*\:(\d+)\]"
        matches = re.findall(pattern, scaffold)  # matches: ['3', '6']
        matches.sort(key=lambda x: int(x))
        old2new = {int(old_idx): new_idx for new_idx, old_idx in enumerate(matches)}
        new2old = {new_idx: int(old_idx) for new_idx, old_idx in enumerate(matches)}
        replacement = "[*:{new_idx}]"
        new_scaffold = re.sub(
            pattern,
            lambda match: replacement.format(new_idx=old2new[int(match.group(1))]),
            scaffold,
        )
        return new_scaffold, new2old, old2new

    @abstractmethod
    def grow(self, *args, **kwargs):
        """
        Abstract method for growing molecules from the scaffold.
        """
        pass

    @abstractmethod
    def clean_tmp_files(self, *args, **kwargs):
        """
        Abstract method for cleaning intermediate files
        """
        pass
