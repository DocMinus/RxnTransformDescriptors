#!/usr/bin/env python
# coding: utf-8
"""test script, creating articifial data, timing and testing functions"""

import time

import pandas as pd

from td_tools.rxntools import clean_smiles_multi, transform_descriptors


class Timer:
    """
    e.g.
    t=Timer()
    codeblock
    t.log()
    """

    # time.perfcounter() better than time.time()?
    def __init__(self):
        self.start = time.time()

    def start(self):
        self.start = time.time()

    def log(self):
        logger = time.time() - self.start
        print("Time log:", logger)

    def milestone(self):
        self.start = time.time()


def main():
    """This test script doesn't use the transform_descriptoer() wrapper, but that isn't necessary for this testing"""
    dataset_size = 10  # change this number to test different timing scenarios
    # we define some faulty/missing compounds, then the output table should have 3 rows less than the input table
    reactant1 = []
    for x in range(0, dataset_size - 1):
        reactant1.append("CCCN")
    reactant1.append("cc")  # incorrect structure
    reactant1.append("CCO")
    reactant1.append("CCCl")
    reactant1.append("CCCl")

    reactant2 = []
    for x in range(0, dataset_size):
        reactant2.append("CCCCO")
    reactant2.append("CC")
    reactant2.append("CCO")
    reactant2.append("CCCl")

    product = []
    for x in range(0, dataset_size):
        product.append("ClCC1=C(B)C(P)=CC(Br)=C1O")
    product.append("cc")  # incorrect structure
    product.append("")  # missing structure
    product.append("CCCl")

    t = Timer()
    g0 = clean_smiles_multi(reactant1)
    g1 = clean_smiles_multi(reactant2)
    g2 = clean_smiles_multi(product)

    TD_numbers = transform_descriptors(g0, g1, g2)

    # create output table
    final_table = pd.DataFrame({"Compound 1": g0, "Compound 2": g1, "Product": g2})
    # filter when empty structures
    final_table = final_table[~((final_table.iloc[:, :3] == "").any(axis=1))]
    final_table = pd.concat([final_table, TD_numbers], axis=1, join="inner")
    # print with some statistics
    print(final_table.info())
    print(final_table.tail())

    print(final_table.shape)
    t.log()  # print the time for descriptor calculation


if __name__ == "__main__":
    main()
