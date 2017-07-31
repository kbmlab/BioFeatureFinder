#!/usr/bin/env python
## Load the packages required

import sys
import os
import argparse
import pandas as pd
import numpy as np
import pybedtools
from pybedtools import BedTool
import warnings
import subprocess
from subprocess import PIPE, Popen
import multiprocessing as mp
from multiprocessing.pool import ThreadPool
from multiprocessing.pool import Pool
import glob

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
pd.options.mode.chained_assignment = None  # default='warn'

##Load the parser for arguments

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        print
        print("The following error ocurred in argument parsing:")
        sys.stderr.write('error: %s\n' % message)
        print
        print(
            "Check the help below. If the error persists, please contact the corresponding author")
        print
        self.print_help()
        sys.exit(2)


## Assign input data as system variables

parser = MyParser(description='')

parser.add_argument('-g', '--gtf', dest="gtf_file",
                    help="GTF file downloaded from Ensembl database. Available at http://www.ensembl.org/info/data/ftp/index.html. Used for analysis of intronic, 3UTR, 5UTR and CDS regions.",
                    required=True)

parser.add_argument('-o', '--outfile', dest="outfile", type=str,
                    help="Name for output files generated.",
                    required=True)

parser.add_argument('-r', '--region', nargs='+', dest="region_list",
                    help="OPTIONAL list of regions for extraction from GTF. Can be passed as a list (Ex. 'three_prime_utr five_prime_utr CDS'), must be an exact match of the 'feature' column in the GTF file. If -r is not used, it detects all features present in the GTF and extract 1 GTF file for each. Default: auto",
                    type=str, default='auto', required=False)

parser.add_argument('--intron', dest="do_introns", action='store_true',
                    help="Use this option to extract intron annotations from the GTF. Requires the presence of genes and exons positions in the annotation. Default: False",
                    required=False, default=False)

parser.add_argument('--split-intron', dest="split_introns", action='store_true',
                    help="Use this option to split intron annotations from the GTF in proximal and distal regions. Requires --intron option. Default: False",
                    required=False, default=False)

parser.add_argument('--analysis', dest="analysis", action='store_true',
                    help="Use this option to perform a preliminary analysis to see in the amount of overlap between a input list of intervals and the regions extracted from the reference GTF. Default: False",
                    required=False, default=False)

parser.add_argument('-i', '--input', dest="input_file", 
                    help="Input list of intervals to overlap with GTF regions. Required for --analysis. Default: False",
                    required=False, default=False)

parser.add_argument("--ncores", dest="ncores", default=(mp.cpu_count() - 1),
                    help="Number of CPU cores used to process feature extraction. Default:(ALL_CORES)-1",
                    type=int, metavar='INT')

args = parser.parse_args()

if args.do_introns is False and args.split_introns is True:
    parser.error("--split-intron requires --intron.")

if args.input_file is False and args.analysis is True:
    parser.error("--analysis requires --input.")
    
def extract_features(feature):
    print
    print "Extracting "+feature+" positions and creating GTF"
    gtf_ref.filter(lambda x: x[2] == feature).saveas('gtf_regions/'+args.outfile+'_'+feature+'.gtf')

print
print "Creating output directory (gtf_regions)"

Popen('mkdir -p gtf_regions', shell=True)
    
print
print "Loading GTF annotation"

gtf_ref = BedTool(args.gtf_file)

if args.region_list == 'auto':
    print
    print "Running in auto mode. Finding region types present in the GTF."
    df = gtf_ref.to_dataframe().dropna()
    feature_list = list(df.feature.value_counts().index)
else:
    feature_list = list(args.region_list)

print
print "Extracting the following regions: "+str(feature_list)
    
p = Pool(args.ncores)
p.map(extract_features, feature_list)

if args.do_introns == True:
    print
    print "Extracting intron positions and generating GTF"
    genes = gtf_ref.filter(lambda x: x[2] == 'gene').saveas()
    exons = gtf_ref.filter(lambda x: x[2] == 'exon').saveas()
    introns = genes.subtract(exons, s=True, nonamecheck=True).saveas()
    introns.saveas('gtf_regions/'+args.outfile+'_intron.gtf')
    
    if args.split_introns == True:
        print
        print "Spliting intron in proximal and distal regions and generating GTF"
        introns_distal = introns.to_dataframe().copy()
        introns_distal.start = introns_distal.start + 500
        introns_distal.end = introns_distal.end - 500
    
        introns_distal_bed = BedTool.from_dataframe(introns_distal
                                                   ).remove_invalid().saveas('gtf_regions/'+args.outfile+'_distal_intron.gtf')
    
        introns_proximal = introns.subtract(introns_distal_bed, s=True, nonamecheck=True
                                            ).saveas('gtf_regions/'+args.outfile+'_proximal_intron.gtf')
    

print
print "GTF files created for each region"

if args.analysis == True:
    input_file = BedTool(args.input_file).sort()
    gtf_list = glob.glob('gtf_regions/*.gtf')
    name_list = [i.split('/', 1)[-1] for i in gtf_list]
    counts = list()
    size = int(input_file.to_dataframe().dropna().shape[0])

    for i in range(len(gtf_list)):
        print
        print "Counting occurences in: "+str(name_list[i])
        ref = BedTool(gtf_list[i])#.sort()
        inter = input_file.intersect(ref, s=True, c=True, nonamecheck=True).to_dataframe()
        counts.append(inter[inter.iloc[:,-1] > 0].shape[0])

    df = pd.DataFrame()
    df['Region'] = name_list
    df['Counts'] = counts
    df['Ratio (counts/input)'] = df.Counts.apply(lambda x: str(round((float(x)/size)*100, 2))+'%')
    df.to_csv('gtf_regions/'+args.outfile+'_gtf_distribution.tsv', sep='\t', index=False)
    
    print
    print "GTF distribution of input intervals"
    print
    print df
else:
    pass
                         
print
print "Zipping output and cleaning temp files"

pybedtools.helpers.cleanup(verbose=False, remove_all=False)
Popen('gzip gtf_regions/*.gtf', shell=True)
                         
print
print "Thank you for using BioFeatureFinder"
print
