#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pytest
import os
from rdkit import Chem
from ta_gen.ocrem.ocrem import mutate_mol
from ta_gen.db import create_db_manager
from pathlib import Path

PKG_DIR = Path(__file__).parents[1]
FRAGMENTATION_EXE = PKG_DIR / "ta_gen" / "bin" / "fragmentation.py"

@pytest.fixture
def test_mol():
    """Create a simple test molecule (toluene)"""
    return Chem.MolFromSmiles('CCO')


@pytest.fixture
def test_smi_file(tmp_path):
    """Create a temporary test.smi file for testing"""
    smi_content = """CC(=O)Oc1ccccc1C(=O)O
CCO
c1ccccc1
"""
    smi_path = tmp_path / "test.smi"
    smi_path.write_text(smi_content)
    return str(smi_path)


@pytest.fixture
def db_manager(test_smi_file, tmp_path):
    """Create a database manager for testing"""
    test_db = tmp_path / "test.db"
    # Create database using fragmentation.py
    if not os.path.exists(test_db):
        os.system(
            f"python {FRAGMENTATION_EXE} --input {test_smi_file} "
            f"--use_db --db_type sqlite --db_path {test_db} --mode 0 --ncpu 10 --radius 3"
        )
    # Create and return database manager
    return create_db_manager('sqlite', db_path=str(test_db))


def test_mutate_mol_basic(test_mol, db_manager):
    """Test basic functionality of mutate_mol"""
    # Test with default parameters
    results = list(mutate_mol(test_mol, db_manager))
    # Should return some results if database has fragments
    assert results == []



