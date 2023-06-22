"""
Collection of variables that might appear at multiple places. 
Easier refactoring, hopefully also better overview(?)
"""
import csv
import os
import re
from rdkit import Chem
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator


PATTERN = "(\[[^\]]+]|Si|Ti|Al|Zn|Pd|Pt|Cu|Br?|Cl?|N|O|S|P|F|I|B|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
REGEX = re.compile(PATTERN)
LOWERCASE_PATTERN = "b|c|n|o|s|p"
REGEX_LOW = re.compile(LOWERCASE_PATTERN)
ELEMENTS = [
    "C",
    "N",
    "O",
    "S",
    "P",
    "F",
    "I",
    "Br",
    "Cl",
    "B",
    "Fe",
    "Te",
    "Se",
    "Sn",
    "Si",
    "Ti",
    "Al",
    "Cu",
    "Zn",
    "Pd",
    "Pt",
]
ELEMENTS_DICT = {key: 0 for key in ELEMENTS}

# specific base Transform descriptors
RDKIT_TD = MolecularDescriptorCalculator(
    [
        "RingCount",
        "NumAliphaticCarbocycles",
        "NumAliphaticHeterocycles",
        "NumAliphaticRings",
        "NumAromaticCarbocycles",
        "NumAromaticHeterocycles",
        "NumAromaticRings",
        "NumRotatableBonds",
        "NumSaturatedCarbocycles",
        "NumSaturatedHeterocycles",
        "NumSaturatedRings",
        "NumSaturatedCarbocycles",
        "NumSaturatedHeterocycles",
        "NumSaturatedRings",
    ]
)

RDKIT_TD_HEADERS = list(RDKIT_TD.GetDescriptorNames())


# fragment based TDs
# read the smarts from a txt based file
TDs_dict = {}
rxnsmarts = os.path.join(os.path.dirname(__file__) + "/rxnsmarts.txt")
with open(rxnsmarts, mode="r") as schmarts:
    reader = csv.reader(schmarts, delimiter="\t")
    TDs_dict = {rows[1]: rows[0] for rows in reader}
TDs_dict.pop("description")
RDKIT_SMARTS = {name: Chem.MolFromSmarts(smarts) for name, smarts in TDs_dict.items()}
