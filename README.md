# BioFeatureFinder: Flexible, unbiased analysis of biological characteristics associated with genomic regions

## Description.

BioFeatureFinder is an algorithm desiged with the goal of uncovering latent biological relationships from sets of genomic coordinates obtained from modern high-throughtput sequencing experiments and subsequent analysis pipelines (i.e. coordinates for protein binding sites or alternatively spliced exons). Designed with flexibility in mind, this algorithm can be used with any species which has a genomic sequence available (gene annotation is highly recommended, but not required) and is compatible with both Python 2.7 and 3.4 in most UNIX-based systems. Due to the modular structure of the algorithm, it can also be easily modified to include novel sources of data and/or analytical functions on demand by reasearchers.

## Scripts:

* ./biofeatures/scripts/analyze_features.py

Main script used for analyzing biological features associated with groups of exons. Uses .bed files of exon coordinates to compare with annotated exons in the data matrix (created by build_datamatrix.py), runs KS statistic to filter non-significant differences and then uses GradientBoost classifier to determine wich features are more “important” in group separation (input vs background).

* ./biofeatures/scripts/build_datamatrix.py

Script used for creating the data matrix for biological features associated with exons and their neighouring regions. Uses a .GTF annotation, genome FASTA, .BW and .BED files as input for BioFeatures. Can be modified to include other functions on demand.

* ./biofeatures/extract_gtf_regions.py

Script used to extract each region (CDS, UTR, Splice Site, Exon, Intron and etc...) from the refference annotation file. By default, it uses the 3rd column of the GTF (Feature) to identify which regions are present in the annotation and extract them. For extracting intronic annotation, it requires the existence of both "gene" and "exon" features in the annotation.

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

### Development install, local changes are reflected in command-line calls (recommended)

    pip install -e .

## Authors

    Felipe E. Ciamponi (GitHub: Fciamponi)
    Michael T. Lovci (GitHub: mlovci)
    Pedro R. S. Cruz
    Katlin B. Massirer

## Funding

BioFeatureFinder was funded by the São Paulo Research Foundation (FAPESP grants 12/00195-3, 13/50724-5, 14/25758-6, 15/25134-5, 16/25521-1).
    
## Testing dataset

For testing of the algorith, we've included a sub-sample of the RBFOX2 RNA-binding protein eCLIP dataset (approx. 10% of the raw dataset). After installing BioFeatureFinder with "pip install -e ." and any dependencies, please also install the following binaries from UCSC utilities directory (http://hgdownload.soe.ucsc.edu/admin/exe/):

    wigToBigWig 
    bigWigMerge 
    bedGraphToBigWig 

Input files for BioFeatureFinder test dataset include the human genome sequence (hg19 assembly) and bigWig files with phastCon scores for conservation. Both of these files are too large to be included in GitHub, so we need to download and process these files to be used by BFF. This can either be done manually (following the steps below) ou using the "get_genome_and_conversation.sh" script located in the "test_data/hg19_data/" folder (WARNING: This step can take a long time. Alternatively, you can download only the genome sequence directly and phastCons scores for multiple alignments of 99 vertebrate genomes (100way) to the human genome already in bigWig format available at: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons.bw).

