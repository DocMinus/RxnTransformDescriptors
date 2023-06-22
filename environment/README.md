# Environment Setup

## Requirements
Python >= 3.9 is required to use the [modern style](https://peps.python.org/pep-0585/) of type annotations.<br>
Recommended: 3.11 (due to increased performance over versions <=3.10)<br>
Modules required are sort of standard for chemistry scripting, rdkit, pandas & numpy, the latter two are nowadays part of a standard conda install.


## Installation with Anaconda/Miniconda
If you nevertheless want a separate environment:<br>
Run the two commands from the root directory.

```shell
conda env create -f ./environment/conda.yaml
conda activate rxn_tds
```

## Installation with Pip
If you already have an environment you want to add this into, then:<br>
Run the command from the root directory

```shell
python -m pip install -r ./environment/requirements.txt
```