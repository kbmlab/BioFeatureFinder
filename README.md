# biofeatures

## Description.


## Warnings


## Scripts:

./biofeatures/scripts/analyze_features.py

Main script used for analyzing biological features associated with groups of exons. Uses .bed files of exon coordinates to compare with annotated exons in the data matrix (created by build_datamatrix.py), runs KS statistic to filter non-significant differences and then uses GradientBoost classifier to determine wich features are more “important” in group separation (input vs background).

./biofeatures/scripts/build_datamatrix.py

Script used for creating the data matrix for biological features associated with exons and their neighouring regions. Uses a .GTF annotation, genome FASTA, .BW and .BED files as input for BioFeatures. 

## Dependency Installation

### With anaconda (recommended):

    conda install -c bioconda pysam pybedtools matplotlib pandas pybedtools scikit-learn matplotlib scipy rpy2
    
### External dependencies (must be available on PATH)
    
    * EMBOSS - http://emboss.sourceforge.net/download/
    * ViennaRNA - https://www.tbi.univie.ac.at/RNA/
    * QGRS Mapper - https://github.com/freezer333/qgrs-cpp
    * BEDTools - http://bedtools.readthedocs.io/en/latest/

## Scripts installation

### Install scripts
    pip install .

### Development install, local changes are reflected in command-line calls

    pip install -e .


## Authors


## Funding


## Pylint

    pylint biofeatures/ > lint_result && git commit -m "ran pylint" lint_result
