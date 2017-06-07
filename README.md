﻿# biofeatures

## Description.


## Warnings


## Scripts:


./biofeatures/scripts/accessory_scripts/convert_annotations.py

./biofeatures/scripts/accessory_scripts/random_sampling_test.py

./biofeatures/scripts/accessory_scripts/reduce_matrix.sh

A script for selecting random lines from the data_matrix.
It's supposed to be used if a computer doesn't have enough RAM
to load the entire matrix as a DataFrame in python.

./biofeatures/scripts/analyze_features.py

./biofeatures/scripts/build_datamatrix.py

Script used for creating the data matrix for biological features associated with exons and their neighouring regions. Uses a .GTF annotation, genome FASTA, .BW and .BED files as input for BioFeatures. Can also use MaxEntScan for calculating splice site score.

## Dependency Installation

### With anaconda (recommended):

    conda install -c bioconda pysam pybedtools matplotlib pandas pybedtools scikit-learn matplotlib scipy rpy2

## Scripts installation

### Install scripts
    pip install .

### Development install, local changes are reflected in command-line calls

    pip install -e .


## Authors


## Funding


## Pylint

    pylint biofeatures/ > lint_result && git commit -m "ran pylint" lint_result
