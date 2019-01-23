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

* ./biofeatures/scripts/analyze_gtf_regions.py

Script used to compare an input BED annotation with multiple GTF references, to indentify preferential regions.

## Dependency Installation

### With anaconda (recommended):

    conda install argparse glob2 matplotlib numpy pandas pybedtools pysam rpy2 scipy seaborn setuptools scikit-learn system
    
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

For testing of the algorith, we've included a sub-sample of the RBFOX2 RNA-binding protein eCLIP dataset (approx. 10% of the raw dataset). For the purpose of this tutorial, we will use only the phastCons scores for multiple alignments of 99 vertebrate genomes to the human genome (100way), which can be obtained directly in bigWig format (for obtaining and processing the 46way dataset, please follow the "Converting multiple wig into single bigWig" section).
    
    #Enter the "hg19_data" inside the "test_data" folder
    cd /full/path/to/BioFeatureFinder/test_data/hg19_data/
    
    #Download the hg19 (GRCh37.p13) from GENCODE ftp website
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz
    
    #Download the 100way phastCons scores from UCSC goldenPath
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons.bw
    
    #(OPTIONAL) Enter the "hg19_annotations" folder and extract 
    #the GTF regions from the Homo sapiens GRCh37.p13 annotation 
    #from Ensembl (with chameleon for UCSC-style)
    cd ./hg19_annotations
    extract_gtf_regions.py -g ./Homo_sapiens.GRCh37.p13.chameleonUCSC.gtf.gz -o grch37 --intron  
    
    #Go back to "test_data" folder
    cd /full/path/to/BioFeatureFinder/test_data/    
    
    #Run region analysis to identify preferential occurence region for input dataset
    analyze_gtf_regions.py \
    -i ./rbfox2_sample.bed \
    -r ./hg19_data/hg19_annotations/grch37.* \
    -l 3pss 3utr 5pss 5utr cds intron \
    -o rbfox2_test_regions
    
    #Run the "build_datamatrix" script to create extract biological features from the files
    build_datamatrix.py \ #load the main script
    -i ./rbfox2_sample.bed \ #set your input set of bed coordinates
    -gen ./hg19_data/GRCh37.p13.genome.fa \ #Load the genomic sequence
    -g ./hg19_data/hg19_annotations/grch37.intron.gtf.gz \ #Load the GTF file 
    -cs ./hg19_data/hg19.100way.phastCons.bw \
    -var ./hg19_data/hg19_var/* \
    --fasta \
    -k 4 5 6 \
    --rnafold \
    --qgrs \
    --keepBED \
    --keepTEMP \
    -o rbfox2.test
    
    #Run the "analyze_features.py" script to identify significant features 
    analyze_features.py \
    -i ./rbfox2_sample.bed \
    -m ./rbfox2.test.datamatrix/rbfox2.test.datamatrix.tsv \
    -o rbfox2.test \
    -p bonferroni \
    -mi 0.01 \
    -c 0.8 \
    -s



## Converting multiple wig into single bigWig

After installing BioFeatureFinder with "pip install -e ." and any dependencies, please also install the following binaries from UCSC utilities directory (http://hgdownload.soe.ucsc.edu/admin/exe/):

    wigToBigWig 
    bigWigMerge 
    bedGraphToBigWig 

Input files for BioFeatureFinder test dataset include the human genome sequence (hg19 assembly) and bigWig files with phastCon scores for conservation. Both of these files are too large to be included in GitHub, so we need to download and process these files to be used by BFF. This can either be done manually (following the steps below) or using some type of automated script.

1 - Download chromosome sizes from UCSC:

    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

2 - Download the md5 checksum for the vertebrate conservation files:

    wget -O md5.vertebrate.txt http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons46way/vertebrate/md5sum.txt

3 - Download each wig file from the md5 and extract them:
    
    cut -f 3 -d " " md5.vertebrate.txt |xargs -P 1 -L 1 -I % wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons46way/vertebrate/% 
    gunzip ./*.wigFix.gz

4 - Remove files without data (can be identified with something like "ls -lah"):
    
    rm chrUn_gl000226.phastCons46way.wigFix

5 - Convert wig files to bigWig:
    
    ls *.wigFix | xargs -P 1 -L 1 -I % wigToBigWig ./% ./hg19.chrom.sizes ./%.bw 

6 - Merge the multiple bigWig files in a single bedGraph and sort it by coordinates:
    
    bigWigMerge ./*.bw hg19.phastCons46way.vertebrates.bg
    bedtools sort -i ./hg19.phastCons46way.vertebrates.bg > ./hg19.phastCons46way.vertebrates.srt.bg
    rm ./*.bw

7 - Convert the bedGraph back into a bigWig binary file and clean up:
    
    bedGraphToBigWig ./hg19.phastCons46way.vertebrates.srt.bg ./hg19.chrom.sizes ./hg19.phastCons46way.vertebrates.bw
    rm ./*.bg

8 - Repeat steps 2 to 7 for other levels of conservation (placentalMammals and primates).
    
## Change logs

* v1.1.4 - analyze_features.py now plots additional visualization options (total of 4: CDF, violin plot, histogram and KDE) in a single multi-page PDF file; analyze_gtf_regions now outputs a clustermap with rudimentary hierarchical clustering using centroid algorithm.
* v1.1.3 - Added early-stopping to GradientBoostClassifier.
* v1.1.2 - Added mutual information score (MI) and linear correlation as feature selection tools.
* v1.1.1 - Fixed issue with hyperparameter tuning when using 10-fold cross validation.
* v1.1.0 - Re-structure of main scripts to support both Python 3.4 and 2.7.
* v1.0.3 - Added test dataset and readme file.
* v1.0.2 - Fixed bug in classifier metrics barchart plotting.
* v1.0.1 - Fixed bug in KS test.
* v1.0.0 - Original BioFeatureFinder algorithm.
