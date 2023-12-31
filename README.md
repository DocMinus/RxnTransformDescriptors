[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![GitHub](https://img.shields.io/github/license/DocMinus/RxnTransformDescriptors)

# Reaction Transform descriptors
Python code to calculate reaction transform descriptors as described in [CHEMRXIV](https://chemrxiv.org/engage/chemrxiv/article-details/649888d41dcbb92a5e8e3475), by [@DocMinus](https://github.com/docminus) and [@DrAlatriste](https://github.com/DrAlatriste). <br>
Not a full fledged package, some scripting know-how necessary to use or incorporate in own code might be necessary.

## Installation
See _environment_ folder.

## Usage
Run the provided script by providing a file with tab/semicolon separated data (also comma or space, though not recommended):<br>
`python 2AB_reaction_TDs.py path/inputfilename`<br>
<br>
You can get help by calling the script using -h: `python 2AB_reaction_TDs.py -h` <br>
<br>
This particular script uses fileformat<br>
_ID reactant1 reactant2 product_<br>
<br>
The script will provide a simple cleaning of the structures; "extreme" broken structures might not get fixed with the provided method.<br>
<br>
Two small test-sets are provided with made up reactions, one of them containing a "faulty" structure to demonstrate correct filtration in the end result. Alternatively, run the _test.py_ script (see below)<br>

## Syntax
If you only want to use the TD function, your script requires the following minimum lines with the smiles as string tuples (even if only a single reaction):
```
    from td_tools.rxntools import transform_descriptors
    
    output_table = transform_descriptors(['smiles_reactant1'],['smiles_reactant2'],['product'])
```
A cleaning function as well as a file reader function is included for larger datasets.<br>
Provided scripts include examples on how to concatenate the structures versus the TDs.<br>
<br>
For quick testing and timing use `Python test.py`.<br>
Not a pytest package, but it nevertheless does the trick for quick demonstrating/testing.<br>
<br>

### Acknowledgments
We would like to thank [@eryl](https://github.com/eryl) for suggestions and help regarding multiprocessing. This allowed processing of large datasets within minutes or even seconds on a standard system, versus previously hours.


