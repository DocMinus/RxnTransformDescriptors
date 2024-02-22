# Environment Setup

## Requirements
Python >= 3.9 is required to use the [modern style](https://peps.python.org/pep-0585/) of type annotations.<br>
Recommended: 3.11 (due to increased performance over vearlier versions)<br>
Modules required are sort of standard for chemistry scripting, rdkit, pandas & numpy, the latter two are nowadays part of a standard conda install. 


## Installation
1. Anaconda/Miniconda
    If you nevertheless want a separate environment:<br>
    Run the two commands from the root directory.

    ```shell
    conda env create -f ./environment/conda.yaml
    conda activate rxn_tds
    ```

    1b. (alternatively) Venv
    Note that venv would also work if you prefer that.

2. Pip
    Now run the requirements with pip into this new environment or into any that you already have.<br>
    Run the command from the root directory

    ```shell
    pip install -r ./environment/requirements.txt
    pip install .
    ```

    The latter installs the rxn_tools into the environment. The example script would work without that, but testing requires that.

## Running Tests
`pytest` is available for testing. See the README.md in /tests.
