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
import matplotlib.pyplot as plt
import seaborn as sns

plt.ioff()

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

if sys.version_info[0] >= 3:
    class BlankLinesHelpFormatter (argparse.HelpFormatter):
        def _split_lines(self, text, width):
            return super()._split_lines(text, width) + ['']
    
    parser = MyParser(description='', 
                      formatter_class=BlankLinesHelpFormatter)

else:
    parser = MyParser(description='')

## Assign input data as system variables

parser.add_argument('-i', '--input', dest="input_list", nargs='+',
                    help="Input list of target genomic regions files for analysis",
                    required=True)

parser.add_argument('-r', '--reference', dest="ref_list", nargs='+',
                    help="Input list of genomic regions reference files for analysis",
                    required=True)

parser.add_argument('-l', '--labels', dest="label_list", nargs='+', type=str,
                    help="Customized labels for region list. Must match the input order of -r list.",
                    required=True)


parser.add_argument('-o', '--outfile', dest="outfile", type=str,
                    help="Name for output files generated.",
                    required=True)


args = parser.parse_args()

if args.label_list is False:
    regions = [i.split('/', -1)[-1] for i in args.ref_list]
else:
    regions = args.label_list
    
name_list = [i.split('/', -1)[-1] for i in args.input_list]

df_list = []

for j in range(len(args.input_list)):
    print
    print("Starting analysis for: "+args.input_list[j])
    input_file  = BedTool(args.input_list[j])#.sort().saveas('sorted.bed')

    counts = list()
    size = int(input_file.to_dataframe().dropna().shape[0])

    for i in range(len(args.ref_list)):
        #print
        #print "Counting occurences in: "+str(regions[i])
        ref = BedTool(args.ref_list[i])#.sort().saveas('sorted2.bed')
        inter = input_file.intersect(ref, f=1, s=True, c=True, nonamecheck=True, sorted=False).to_dataframe()
        counts.append(inter[inter.iloc[:,-1] > 0].shape[0])

    df = pd.DataFrame()
    df['Region'] = regions
    df[name_list[j]] = counts
    df.set_index('Region', inplace=True)
    df_list.append(df)

os.makedirs('./'+args.outfile+'_region_analysis/', 
            exist_ok=True)
    
df_cat = pd.concat(df_list, axis=1).T
df_cat.to_excel('./'+args.outfile+'_region_analysis/'+args.outfile+'.xlsx')
df_cat.to_csv('./'+args.outfile+'_region_analysis/'+args.outfile+'.csv')

if len(args.input_list) >= 2:
    
    print("Plotting clustermap")
    print()
    
    g = sns.clustermap(df_cat, 
                       z_score=0,
                       metric='seuclidean',
                       method='single',
                       linewidths=.5,
                       annot=True, fmt=".2f",
                       cmap="coolwarm",
                       col_cluster=True,
                       row_cluster=True,
                       figsize=(np.array(df_cat.T.shape)/1.5))
    
    g.ax_heatmap.set_title('Clustermap of Regions x Input', 
                           size=20, position=(0.5,1.3))
    
    g.ax_heatmap.set_xlabel('Regions', size=20)
    g.ax_heatmap.set_ylabel('Input', size=20)
    g.ax_heatmap.tick_params(axis='both', which='major', labelsize=14)
    
    plt.savefig('./'+args.outfile+'_region_analysis/clustermap.pdf',
                dpi=300, bbox_inches='tight')
    
    plt.savefig('./'+args.outfile+'_region_analysis/clustermap.svg',
                dpi=300, bbox_inches='tight')
    
    plt.savefig('./'+args.outfile+'_region_analysis/clustermap.jpg',
                dpi=300, bbox_inches='tight', quality=95)
    
    plt.close()
else:
    pass

pybedtools.helpers.cleanup(verbose=False, remove_all=False)

print("Finished region analysis. Thank you for using BioFeatureFinder")
print()
