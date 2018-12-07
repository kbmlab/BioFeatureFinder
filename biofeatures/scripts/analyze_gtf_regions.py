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
import scipy
from scipy.cluster.hierarchy import set_link_color_palette
from itertools import compress, product

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

parser.add_argument('-f', '--fraction', dest="f_val", type=float,
                    help="Fraction of the input that must overlap with reference region (must be float from 0 to 1). Default: False",
                    required=False, default=False)


parser.add_argument('-o', '--outfile', dest="outfile", type=str,
                    help="Name for output files generated.",
                    required=True)


args = parser.parse_args()

def combinations(items):
    return ( set(compress(items,mask)) for mask in product(*[[0,1]]*len(items)) )

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
        inter = input_file.intersect(ref, f=args.f_val, s=True, c=True, nonamecheck=True, sorted=False).to_dataframe()
        counts.append(inter[inter.iloc[:,-1] > 0].shape[0])

    df = pd.DataFrame()
    df['Region'] = regions
    df[name_list[j]] = counts
    df.set_index('Region', inplace=True)
    df_list.append(df)

os.makedirs('./'+args.outfile+'_region_analysis/', 
            exist_ok=True)
    
df_cat = pd.concat(df_list, axis=1).T

if len(args.input_list) >= 2:
    
    print("Plotting clustermap")
    print()
        
    g = sns.clustermap(df_cat, 
                       z_score=0,
                       #standard_scale=0,
                       metric='euclidean',
                       method='centroid')
    
    # plt.savefig('./' + args.prefix + '.analysis/classifier_metrics/mean_relative_importance.pdf',
    #             dpi=300, bbox_inches='tight')
    
    plt.close()
    
    n_colors = len(list(combinations(range(df_cat.shape[1]))))
    
    set_link_color_palette(list(sns.color_palette("cubehelix", n_colors).as_hex()))
    
    den = scipy.cluster.hierarchy.dendrogram(g.dendrogram_row.linkage,
                                             labels = df_cat.index,
                                             color_threshold=(0.7*max(g.dendrogram_row.linkage[:,2])),
                                             )
    
    n_clusters = len(set(list(filter(lambda k: '#' in k, den['color_list']))))
    
    set_link_color_palette(list(sns.color_palette("cubehelix", n_clusters).as_hex()))
    
    den = scipy.cluster.hierarchy.dendrogram(g.dendrogram_row.linkage,
                                             labels = df_cat.index,
                                             color_threshold=(0.7*max(g.dendrogram_row.linkage[:,2])),
                                             )
    
    plt.close()
    
    from collections import defaultdict
    
    def get_cluster_classes(den, label='ivl'):
        cluster_idxs = defaultdict(list)
        for c, pi in zip(den['color_list'], den['icoord']):
            for leg in pi[1:3]:
                i = (leg - 5.0) / 10.0
                if abs(i - int(i)) < 1e-5:
                    cluster_idxs[c].append(int(i))
    
        cluster_classes = {}
        for c, l in cluster_idxs.items():
            i_l = [den[label][i] for i in l]
            cluster_classes[c] = i_l
    
        return cluster_classes
    
    clusters = get_cluster_classes(den)
    
    cluster = []
    for i in df_cat.index:
        included=False
        for j in clusters.keys():
            if i in clusters[j]:
                cluster.append(j)
                included=True
        if not included:
            cluster.append(None)
            
    df_cat['cluster'] = cluster
    
    g = sns.clustermap(df_cat.drop("cluster",1), 
                       z_score=0,
                       #standard_scale=0,
                       metric='euclidean',
                       method='centroid',
                       linewidths=.5,
                       annot=True, fmt=".2f",
                       cmap="coolwarm",
                       col_cluster=True,
                       row_cluster=True,
                       figsize=(np.array(df_cat.T.shape)/1.5),
                       row_colors=df_cat.cluster)
    
    g.ax_heatmap.set_title('Clustermap of Input x Region', 
                           size=20, position=(0.5,1.3))
    
    g.ax_heatmap.set_xlabel('Regions', size=20)
    g.ax_heatmap.set_ylabel('Input', size=20)
    g.ax_heatmap.tick_params(axis='both', which='major', labelsize=14)
    g.ax_row_colors.tick_params(axis='both', which='major', labelsize=14)
    
    xlabels = g.ax_heatmap.get_xticklabels()
    g.ax_heatmap.set_xticklabels(xlabels, rotation=90)

    ylabels = g.ax_heatmap.get_yticklabels()
    g.ax_heatmap.set_yticklabels(ylabels, rotation=0)
    
    g.savefig('./'+args.outfile+'_region_analysis/clustermap.pdf',
                dpi=300, bbox_inches='tight')
    
    g.savefig('./'+args.outfile+'_region_analysis/clustermap.svg',
                dpi=300, bbox_inches='tight')
    
    g.savefig('./'+args.outfile+'_region_analysis/clustermap.jpg',
                dpi=300, bbox_inches='tight', quality=95)
    
    plt.close()
else:
    pass


df_cat.to_excel('./'+args.outfile+'_region_analysis/'+args.outfile+'.xlsx')
df_cat.to_csv('./'+args.outfile+'_region_analysis/'+args.outfile+'.csv')

pybedtools.helpers.cleanup(verbose=False, remove_all=False)

print("Finished region analysis. Thank you for using BioFeatureFinder")
print()
