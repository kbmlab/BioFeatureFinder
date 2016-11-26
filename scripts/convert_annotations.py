## Load the packages required

import sys
import os
import argparse
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import pybedtools
from pybedtools import BedTool
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 
pd.options.mode.chained_assignment = None  # default='warn'

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        print
        print "The following error ocurred in argument parsing:"
        sys.stderr.write('error: %s\n' % message)
        print
        print "Check the help below and try to fix the arguments. If the error persists, please contact the corresponding author"
        print
        self.print_help()
        sys.exit(2)

## Assign input data as system variables

parser = MyParser(description='')

parser.add_argument("--ucsc-to-ensembl",dest="ucsc_to_ensembl", 
                  action="store_true", default=False,
                  help="Makes files obtained in UCSCs database compatible with Ensembl chr annotation.")

parser.add_argument("--ensembl-to-ucsc",dest="ensembl_to_ucsc", 
                  action="store_true", default=False,
                  help="Makes files obtained in Ensembl database compatible with UCSCs chr annotation.")

parser.add_argument('-f', '--file', dest="conv_file", 
                  help="File to be converted", metavar="FILE", required=True)

parser.add_argument('-o', '--output', dest="out", 
                  help="Name of output file", metavar="OUT", required=True)

args = parser.parse_args()

if args.ensembl_to_ucsc == True:
    ##Convert Ensembl to UCSC
    ens = pd.concat(BedTool(args.conv_file).to_dataframe(iterator=True, chunksize=10000), ignore_index=True).dropna()
    ens = ens[~ens[ens.columns[0]].str.contains('GL|KI').fillna(False)]
    ens[ens.columns[0]] = 'chr'+ens[ens.columns[0]].astype(str)
    ens[ens.columns[0]] = ens[ens.columns[0]].str.replace('chrMT','chrM')
    columns = ens.select_dtypes(include=['float']).columns.tolist()
    ens[columns] = ens[columns].astype(int)
    ens.to_csv(args.out, compression='gzip', sep='\t', index=False, header=False)
    
elif args.ucsc_to_ensembl == True:
    ##Convert UCSC to Ensembl
    ucsc = pd.concat(BedTool(args.conv_file).to_dataframe(iterator=True, chunksize=10000), ignore_index=True).dropna()
    ucsc[ucsc.columns[0]] = ucsc[ucsc.columns[0]].str.replace('chrM','chrMT')
    ucsc[ucsc.columns[0]] = ucsc[ucsc.columns[0]].apply(lambda x: x.split('chr')[1])
    columns = ucsc.select_dtypes(include=['float64']).columns.tolist()
    ucsc[columns] = ucsc[columns].astype(int)
    ucsc.to_csv(args.out, compression='gzip', sep='\t', index=False, header=False)
    
else:
    print "No conversion option selected. Please use --ucsc-to-ensembl or --ensembl-to-ucsc."