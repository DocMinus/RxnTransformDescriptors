#!/usr/bin/env python
# coding: utf-8
""" test script, creating articifial data and testing the TDs calculation for reactions
    only tests the combination and final outcome, not the individual functions
    2024-02-22; DocMinus
"""

import pandas as pd
import pytest

from td_tools.rxntools import clean_smiles_multi, transform_descriptors


def test_clean_smiles_multi_and_transform_descriptors():
    dataset_size = 4  # number of compounds
    # we define some faulty/missing compounds, then the output table should have 3 rows less than the input table
    reactant1 = ["CCCN" for _ in range(dataset_size - 1)]
    reactant1.append("cc")  # incorrect structure
    reactant1.append("CCO")
    reactant1.append("CCCl")
    reactant1.append("CCCl")
    total_dataset_size = len(reactant1)

    reactant2 = ["CCCCO" for _ in range(dataset_size)]
    reactant2.append("CC")
    reactant2.append("CCO")
    reactant2.append("CCBr")

    product = ["ClCC1=C(B)C(P)=CC(Br)=C1O" for _ in range(dataset_size)]
    product.append("cc")  # incorrect structure
    product.append("")  # missing structure
    product.append("CCI")

    """ a total of 3 rows faulty rows, from bottom of created table it would 2nd, 3rd and 4th last row."""

    g0 = clean_smiles_multi(reactant1)
    g1 = clean_smiles_multi(reactant2)
    g2 = clean_smiles_multi(product)

    TD_numbers = transform_descriptors(g0, g1, g2)
    print(f"{TD_numbers.shape = }, {TD_numbers.shape[1] = }")

    final_table = pd.DataFrame({"Compound 1": g0, "Compound 2": g1, "Product": g2})
    final_table = final_table[~((final_table.iloc[:, :3] == "").any(axis=1))]
    final_table = pd.concat([final_table, TD_numbers], axis=1, join="inner")

    # check if the final table is as expected (3 rows less than the input table)
    assert (
        final_table.shape[0] == total_dataset_size - 3
    ), "Number of rows in the final table is not as expected."

    # Check that the 2nd, 3rd, and 4th last rows have been removed
    removed_indices = [
        total_dataset_size - 2,
        total_dataset_size - 3,
        total_dataset_size - 4,
    ]
    for index in removed_indices:
        assert (
            index not in final_table.index
        ), f"Row {index} should have been removed but is still in the final table."
