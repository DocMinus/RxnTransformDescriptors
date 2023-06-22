[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![GitHub](https://img.shields.io/github/license/docminus/RxnTransformDescriptors)

# Reaction Transform descriptors
Python code to calculate reaction transform descriptors as described in CHEMRXIV **update when available**

## Installation
See _environment_ folder.

## Usage
Run the provided script by providing your tabe/semicolon (also comma or space, though not recommended):<br>
`python 2AB_reaction_TDs.py path/inputfilename`<br>
You can get help by calling the script using -h: `python 2AB_reaction_TDs.py -h` <br>
This particular script uses fileformat<br>
_ID reactant1 reactant2 product_<br>
The script will provide a simple cleaning of the structures; "extreme" broken structures might not get fixed with the provided method.<br>
<br>
Two small test-sets are provided with made up reactions, one of them containing a "faulty" structure to demonstrate correct filtration in the end result.

### Syntax
If you only want to use the TD function, your script requires the following minimum lines with the smiles in tuples (even if only a single one):
```
    from td_tools.rxntools import transform_descriptors
    
    output_table = transform_descriptors(['smiles_reactant1'],['smiles_reactant2'],['product'])
```
A cleaning function as well as a file reader function is included for larger datasets.<br>
Provided scripts included examples on how to concatenate the structures versus the numbers.<br>
<br>
For quick testing and timing use `test.py` script.<br>
Not a pytest package but does the trick.<br>
<br>

