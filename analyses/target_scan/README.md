# target_scan

Predict targets of the two miRNAs using TargetScan

## How to use

* Use the `results` directory during exploration but remember to save milestone files in the project directory.
* Exclude large data from version control by adding it to `.gitignore`

## Custom Environments

If you create a new conda environment for this analysis. Please remember save
your environment to a file after each change:

`mamba env export > environment.yml`

You can install the python package in this directory to your local conda environment with:

`pip install -e .`

You will need to manually add large data files to .gitignore to prevent it from syncing to
version control.

## Authors

* felixl