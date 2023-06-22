#!/usr/bin/env python
# coding: utf-8
"""
V2.1.3 (Mar. 08, 08:00:00 2023)
Update: 2023-06-22 (cleanup for ChemRxiv submission)

@author: Alexander Minidis (DocMinus)
Purpose: TDs from csv
"""

import argparse
import os

import pandas as pd

from td_tools.rxntools import clean_smiles_multi, read_rct2pd, transform_descriptors


def main():
    #
    # Initialization
    #
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            """\
        Calculate transform and elemental analysis descriptors for reactions.
        Basic cleaning included.
        Format of file: Smiles with columns like this:
            ID  Component1 Component2 Product
        preferably tab or semicolon separated. 
        """
        ),
    )
    parser.add_argument("file", type=str)
    args = parser.parse_args()
    if not ("data_path" in locals()):
        file_in = os.path.basename(args.file)
        data_path = os.path.dirname(args.file)
        if data_path == "" or data_path == ".":
            data_path = os.getcwd()

    file_base = os.path.splitext(file_in)[0]
    smiles_input = os.path.join(data_path + "/" + file_in)
    final_output_tsv = os.path.join(data_path + "/" + file_base + "_TD.tsv")
    final_output_pkl = os.path.join(data_path + "/" + file_base + "_TD.pkl")

    #############################################################################
    # Read
    print("\nReading molecules from: ", smiles_input)
    in_rct_df = read_rct2pd(smiles_input)

    # Clean
    cmpd1_smi = clean_smiles_multi(in_rct_df["compound1"].to_list())
    cmpd2_smi = clean_smiles_multi(in_rct_df["compound2"].to_list())
    prod_smi = clean_smiles_multi(in_rct_df["product"].to_list())

    # Calculate TDs
    transforms_descriptors = transform_descriptors(cmpd1_smi, cmpd2_smi, prod_smi)

    # for output create table with structures and combine with calculated TDs
    _df = pd.DataFrame(
        {"Compound 1": cmpd1_smi, "Compound 2": cmpd2_smi, "Product": prod_smi}
    )
    # In addition: filter when empty structures
    _df = _df[~((_df.iloc[:, :3] == "").any(axis=1))]
    # The three tables are concatenated to one
    final_table = pd.concat(
        [in_rct_df["ID"], _df, transforms_descriptors], axis=1, join="inner"
    )
    # output (optional)
    print(final_table.tail())
    #############################################################################
    # Write pickle & csv file
    print("\nWriting to file: ", final_output_pkl)
    final_table.to_pickle(final_output_pkl)
    print("\nWriting to file: ", final_output_tsv)
    final_table.to_csv(
        final_output_tsv, sep="\t", header=True, index=False, quoting=None, decimal="."
    )


if __name__ == "__main__":
    main()
