"""
Initial Version: Nov, 2022
Version 2.1.3
Update: 2023-06-22
@author: Alexander Minidis (DocMinus)

Copyright (c) 2022-2023 DocMinus
Changelog: 
V1 -> V2: multiprocessing with help by @eryl (github)
Finalization w. further optimization and removal of V1 code. Black-ened
Also added csv reader to pd (copy from chemtools)
Added transfrom_descriptors() as wrapper for all descriptor functions, reducing number of lines in main.py
"""

import multiprocessing
from collections import defaultdict

# general imports
import pandas as pd

# RDkit stuff
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

# my imports
from .ids import (
    ELEMENTS_DICT,
    RDKIT_SMARTS,
    RDKIT_TD,
    RDKIT_TD_HEADERS,
    REGEX,
    REGEX_LOW,
)

# Helpful, else way too many RDkit details in output
RDLogger.logger().setLevel(RDLogger.CRITICAL)


def is_smiles(text: str) -> bool:
    """
    Determine if input possibly is a smiles or just a regular text
    Tokenize the input and compare the input length vs token length
    If not the same, then a regular word or phrase.
    based on: https://github.com/pschwllr/MolecularTransformer
        Input:
            string
        Output:
            boolean, True if smiles, False if regular text
    """
    # This simple pre-check seems to work for catching plain numbers as ID (else detected as smiles)
    # normally isnumeric could be used, but if it is a number, it's from numpy and throws error.
    if not isinstance(text, str):
        return False

    tokens = [token for token in REGEX.findall(text)]
    length = len("".join(tokens))
    # empty could give len 0!
    # >= probably not necessary, but safer than ==
    if length and length >= len(text):
        return True
    else:
        return False


def read_rct2pd(file_name: str, delimiter="\t", header=None) -> pd.DataFrame:
    """
    Assumes the file to be structured according to:

    ID in column zero
    SMILES for A+B->C in columns one, two and three. Any additional columns are dropped.

    :param file_name: (incl path)
    :param delimiter:  optional
    :param header: not really used due to "auto" detection
    :return: pandas table w ID in first column & reaction components in 2,3 & 4
    """

    with open(file_name) as file_in:
        first_line = file_in.readline().split()
    if len(first_line) < 4:
        print("---File Error: minimum 4 columns: 'ID and A + B & Product' ---")
        return

    # column 1 as rep for compound or not
    if is_smiles(first_line[1]):
        df = pd.read_csv(file_name, sep=delimiter, header=None, dtype=None)
        df.rename(
            columns={
                0: "ID",
                1: "component1",
                2: "component2",
                3: "product",
            },
            inplace=True,
        )
    else:
        df = pd.read_csv(file_name, sep=delimiter, header=0, dtype=None)

    # drop columns that aren't ID or part of reaction i.e. keep only first 4 columns
    if len(df.columns) > 4:
        df.drop(df.columns[4:], axis=1, inplace=True)

    # coerce ID to string in case it still is handled as number. to prevent later errors.
    df.iloc[0:1] = df.iloc[0:1].astype("str")
    return df


def rd_clean(workpackage: str) -> list:
    """RDKIT minimal cleaning that takes care of majority of cases.\n
    Will not clean extremely 'dirty' sources.\n
    Incorrect smiles will be returned as empty.\n
    Returns input ID & smiles, not object, better for memory/multiprocessing.

    in: single smiles string
    out: id & cleaned smiles string
    """
    _id, smiles = workpackage
    mol = Chem.MolFromSmiles(smiles, sanitize=False)

    try:
        Chem.SanitizeMol(mol)
    except ValueError as _e:
        print(f"ID: {_id}_{_e}")
        return _id, ""

    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(
        mol,
        sanitizeOps=(
            Chem.SANITIZE_ALL ^ Chem.SANITIZE_CLEANUP ^ Chem.SANITIZE_PROPERTIES
        ),
    )
    # canonicalize:
    mol = rdMolStandardize.Normalize(mol)
    return _id, Chem.MolToSmiles(mol)


def clean_smiles_multi(schmiles: list) -> tuple:
    """parallel cleanining of a list of smiles.\n
    Includes an ID, not really necessary, but helps in later analysis downstream.

    in: list of smiles
    out: list of cleaned smiles (for now, without ID)
    """
    print("...Cleaning compounds...")
    work = [(i, smiles) for i, smiles in enumerate(schmiles)]
    with multiprocessing.Pool() as pool:
        processed_smiles = pool.map(rd_clean, work)

    # tuple or list, doesn't really matter(?)
    return tuple(smi[1] for smi in processed_smiles)


def element_count(smiles: str) -> dict:
    """
    Creates a "sum formula" as dictionary, no hydrogens though.\n
    Doesn't consider salts, ignores bondtypes. But: doesn't need RDKIT.
        Input:
            string containing smiles
        Output:
            dictionary containing sumformula, here: sorted, including 0 atoms, useful for EA descriptors
    """

    elements = [token for token in REGEX.findall(smiles)]
    for i in range(len(elements)):
        if REGEX_LOW.findall(elements[i]):
            elements[i] = elements[i].upper()

    element_count = [elements.count(ele) for ele in elements]
    formula = dict(zip(elements, element_count))
    # final_dict = {"smiles": smiles}
    # final_dict.update(ELEMENTS_DICT)
    # return {key: formula.get(key, final_dict[key]) for key in final_dict}
    return {key: formula.get(key, ELEMENTS_DICT[key]) for key in ELEMENTS_DICT}


def elemental_tds_multi(schmiles: list) -> pd.DataFrame:
    """
    Calculates the Elemental Analysis based TDs for a reaction.\n

    :param _list: pandas with the (cleaned) input reaction
    :return: pandas with the TDs based on elemental analysis
    """
    print("...Calculation EA based descriptors...")

    with multiprocessing.Pool() as pool:
        processed_smiles = pool.map(element_count, schmiles)

    # return processed_smiles
    return pd.DataFrame(data=processed_smiles).astype("int16")


def properties_calc(schmiles: str) -> list:
    """
    RDKIT (selected) descriptor calculations \n
    Best to do cleaning/error checking before since not included.
        Input:
            smiles string
        Output:
            dictionary containing single mol based TDs
    """

    return RDKIT_TD.CalcDescriptors(Chem.MolFromSmiles(schmiles))


def rdkit_descriptors_multi(schmiles: list) -> pd.DataFrame:
    """
    calculates RDkit descriptors.\n
    Best to do error/cleaning of smiles before using this.

    :param list of smiles (! no cleaning/check for error)
    :return: pandas df containing calculated properties (no structures, but now with ID)
    """
    print("...Calculating RDKit descriptors...")
    with multiprocessing.Pool() as pool:
        processed_smiles = pool.map(properties_calc, schmiles)

    named_properties = dict(zip(RDKIT_TD_HEADERS, zip(*processed_smiles)))
    return pd.DataFrame(data=named_properties).astype("int16")


def smarts_count(schmiles: str) -> dict:
    """
    Counts fragment based transforms.\n
    No cleaning/error checking of smiles\n
    Smarts stem from text file rxnsmarts.txt.
        Input:
            smiles string
        Output:
            dictionary containing single mol based TDs
    """

    rdmol = Chem.MolFromSmiles(schmiles)
    matches = {}
    for name, smarts in RDKIT_SMARTS.items():
        matches[name] = len(rdmol.GetSubstructMatches(smarts))

    return matches


def smarts_descriptors_multi(schmiles: list) -> pd.DataFrame:
    """
    Parallel optimized smarts calculation. Calls smarts_count.
        Input:
            list of smiles (not rd objects)
        Output:
            results in a DataFrame.
    """
    print("...Calculating smarts based descriptors...")

    with multiprocessing.Pool() as pool:
        processed_smiles = pool.map(smarts_count, schmiles)
        # processed_smiles = list(pool.imap_unordered(match_smarts, mollist, chunksize=chunksize))

    columns = defaultdict(list)
    for smi in schmiles:
        columns["smiles"].append(smi)

    for smart_matches in processed_smiles:
        for name, match in smart_matches.items():
            columns[name].append(match)

    return pd.DataFrame(columns).iloc[:, 1:].astype("int16")


def table_delta(
    cmpd1: pd.DataFrame, cmpd2: pd.DataFrame, prod: pd.DataFrame
) -> pd.DataFrame:
    """very simple pandas table sub/add. Input order important!"""
    return prod.sub(cmpd1.add(cmpd2))


def transform_descriptors(
    cmpd1_smi: tuple, cmpd2_smi: tuple, prod_smi: tuple
) -> pd.DataFrame:
    """Wrapper calling all necessary functions to calculate TDs"""

    eas0 = elemental_tds_multi(cmpd1_smi)
    eas1 = elemental_tds_multi(cmpd2_smi)
    eas2 = elemental_tds_multi(prod_smi)
    ea_final = table_delta(eas0, eas1, eas2)

    rdk0 = rdkit_descriptors_multi(cmpd1_smi)
    rdk1 = rdkit_descriptors_multi(cmpd2_smi)
    rdk2 = rdkit_descriptors_multi(prod_smi)
    rdk_final = table_delta(rdk0, rdk1, rdk2)

    smarts1 = smarts_descriptors_multi(cmpd1_smi)
    smarts2 = smarts_descriptors_multi(cmpd2_smi)
    smarts3 = smarts_descriptors_multi(prod_smi)
    smarts_final = table_delta(smarts1, smarts2, smarts3)

    # The three tables are concatenated to one, the TD table
    return pd.concat([ea_final, rdk_final, smarts_final], axis=1, join="inner")
