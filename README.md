[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![GitHub](https://img.shields.io/github/license/DocMinus/RxnTransformDescriptors)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/DocMinus/RxnTransformDescriptors)

# Reaction Transform descriptors
Python code to calculate reaction transform descriptors as described in [CHEMRXIV](https://chemrxiv.org/engage/chemrxiv/article-details/649888d41dcbb92a5e8e3475), by [@DocMinus](https://github.com/docminus) and [@DrAlatriste](https://github.com/DrAlatriste). <br>

## Installation
See _environment_ folder.
Updated the installation with a setup file to enable the tools to be part of ones Python environment. Testing has also been added. 

## Example Usage
Run the example script by providing a file with tab/semicolon separated data (also comma or space, though not recommended):
```shell
python AB2C_reaction_TDs_example.py inputfilename
```
<br>
You can get help by calling the script using -h: `python AB2C_reaction_TDs_example.py -h` <br>
<br>
This particular script expects the input order of the file as<br>

_ID reactant1 reactant2 product_ <br>
<br>
Simple cleaning of structures is included; "extreme" broken structures might not get fixed with the provided method.
<br>
Two small test-sets are provided with made up reactions, one of them containing a "faulty" structure to demonstrate correct filtration in the output result. <br>Execute via:  `python AB2C_reaction_TDs_examples.py ./datsets/testreactions.tsv`<br>

## Syntax
If you only want to use the TD function, your script requires the following minimum lines with the smiles as string tuples (even if only a single reaction):
```shell
from td_tools.rxntools import transform_descriptors

output_table = transform_descriptors(['smiles_reactant1'],['smiles_reactant2'],['product'])
```
A cleaning function as well as a file reader function is included for larger datasets:
```shell
from td_tools.rxntools import clean_smiles_multi, read_rct2pd  
```
The provided script includes examples on how to concatenate the structures versus the TDs.<br>

## Testing
Python testing has been added instead of the previous test.py, see the README.md under /tests.<br>
<br>

### Acknowledgments
We would like to thank [@eryl](https://github.com/eryl) for suggestions and help regarding multiprocessing in the original build. This allowed processing of large datasets within minutes or even seconds on a standard system, versus previously hours.<br>
Currently this has been changed to joblib instead, which seems a bit more stable and faster in this particular context.

### Updates
* setup.py for install as package
* testing added
* switch from multiparallel to joblib
* releases introduced; version number reflects version number of tool.
