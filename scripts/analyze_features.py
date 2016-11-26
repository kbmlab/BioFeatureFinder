#### Load the required packages

import sys
import argparse
from argparse import ArgumentParser
import os
import time
import pandas as pd
import numpy as np
import pybedtools
from pybedtools import BedTool
import scipy
from scipy import stats
import readline
import rpy2
import rpy2.robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import itertools
from itertools import cycle
import matplotlib
import matplotlib.pyplot as plt
import sklearn
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble.partial_dependence import plot_partial_dependence
from sklearn.ensemble.partial_dependence import partial_dependence
from sklearn.metrics import adjusted_mutual_info_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import mean_squared_error
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
%matplotlib inline
pd.options.mode.chained_assignment = None  # default='warn'


##Load the parser for arguments

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

parser.add_argument('-b', '--bed', dest="bed_file", 
                  help="BED file with exons/regions of interest.", required=True)

parser.add_argument('-ref', '--refference_bed', dest="refference_bed", 
                  help="BED file with annotated exons created by 'build_datamatrix.py'", required=True)

parser.add_argument('-m', '--matrix', dest="matrix", 
                  help="Data matrix with biological features created by 'build_datamatrix.py'", required=True)

parser.add_argument('-o', '--outfile', dest="prefix", 
                  help="prefix for use on the output files", metavar="prefix", required=True)

parser.add_argument("-padj", '--p_adjust', dest="padj", default='bonferroni', help="Type of p-value correction used after Kolmogorov-Smirnov test, available options are: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' and 'none'. Default:'bonferroni'", type=str, metavar='padj', required=False)

parser.add_argument('-filter', '--filter_columns', dest="filter_out", default=False,
                  help="Text file containing a comma-separated list with names of the columns to be removed from the dataframe in the analysis. Default: False", metavar="filter_out.txt", required=False)

parser.add_argument('-select', '--select_columns', dest="filter_in", default=False,
                  help="Text file containing a comma-separated list with names of the columns in dataframe to be used in the analysis. Default: False", metavar="filter_in.txt", required=False)

parser.add_argument("--ncores", dest="ncores", default=(mp.cpu_count()-1), help="Number of CPU cores used to process downtream/upstream exon search. Default:(ALL_CORES)-1", type=int, metavar='INT')

parser.add_argument("--debug",dest="debug", metavar='INT',
                    type=int, default=False, 
                    help="Only seaches for the first N transcripts. Useful for debugin and checking if the code is working properly. Default: False")

args = parser.parse_args()

##Define the functions which will be used during the analysis

def group_matrices_one_sample(bt, bt_a, matrix, F=1, f=1, s=True):
    
    int_a = bt.intersect(bt_a, 
                         s=s, 
                         F=F, 
                         f=f, 
                         sorted=True).to_dataframe().drop_duplicates()
    
    matrix['group'] = 0
    matrix['group'][matrix['name'].isin(int_a['name'])] = 1
    
    return matrix


def get_statistical_data_for_features(df, correction_type=args.padj):
    features = pd.DataFrame()
    features['Feature'] = df.drop('group',1).columns
    
    pairwise = list(itertools.combinations(df.group.drop_duplicates(), 2))
    
    statsR = importr('stats')
    
    for i in range(len(pairwise)):
        print "Calculating Kolmogorov-Smirnov test for: Group "+str(pairwise[i][0])+" vs Group "+str(pairwise[i][1])
        print
        
        features['ks_'+str(pairwise[i][0])+'_vs_'+str(pairwise[i][1])] = features['Feature'].apply(
            lambda x: stats.ks_2samp(df[df['group'] == pairwise[i][0]][x].astype(float),
                                     df[df['group'] == pairwise[i][1]][x].astype(float)))
    
        features['pval_'+str(pairwise[i][0])+'_vs_'+str(pairwise[i][1])
                ] = features['ks_'+str(pairwise[i][0])+'_vs_'+str(pairwise[i][1])].apply(lambda x: x[:][1])
        
        features['ks_'+str(pairwise[i][0])+'_vs_'+str(pairwise[i][1])
                ] = features['ks_'+str(pairwise[i][0])+'_vs_'+str(pairwise[i][1])].apply(lambda x: x[:][0])
        
        print "Adjusting P-value scores with "+str(correction_type)+" method"
        print
        
        features['adj_pval_'+str(pairwise[i][0])+'_vs_'+str(pairwise[i][1])
                ] = statsR.p_adjust(FloatVector(features['pval_'+str(pairwise[i][0])+'_vs_'+str(pairwise[i][1])]), 
                                   method = correction_type)
                                       
        print "Done with KS test for: Group "+str(pairwise[i][0])+" vs Group "+str(pairwise[i][1])
        print
    return features

def cdf(data, bins=50):
    #data = np.ma.masked_array(data, np.isnan(data))
    data = data.astype(float)
    minimum = np.min(data) - .000001
    maximum = np.max(data) + .000001
    pos = np.linspace(minimum, maximum, bins + 1)
    xs = np.linspace(minimum, maximum, bins + 1)[:-1]
    ys = np.linspace(minimum, maximum, bins + 1)[1:]
    ecdf = np.ndarray(shape=(bins + 1, 1))
    ecdf[0] = 0
    cumSum = 0
    for i, (x, y) in enumerate(zip(xs, ys)):
        region = len(data[np.where((data >= x) & (data < y))])
        cumSum += region / float(len(data))
        ecdf[i + 1] = cumSum
    return pos, ecdf

def plot_cdf(data, bins=50, ax=None, **plotting_args):
    if ax is None:
        ax = plt.gca()
    x, y = cdf(data, bins=bins)
    ax.plot(x,y, **plotting_args)
    
##Load the datamatrix generatade by "buildadatamatrix.py"

print
print "Loading datamatrix"
print

matrix = pd.concat(pd.read_table(args.matrix,iterator=True, chunksize=10000), 
                          ignore_index=True).set_index('name').drop_duplicates().reset_index()

##Load the bed file created with all exons

print "Loading annotated exons in refference file"
print 

all_exons = BedTool(args.refference_bed).sort()

##Load the dataset of interest to be analyzed (needs to be a SORTED BED)

print "Load selected exons for biological feature analysis"
print

sel_exons = BedTool(args.bed_file).sort()

##Intersect the exons found in the analysis to get groups 1 (positive) and 0 (negative) in the matrix

print "Finding input exons in the matrix and selecting groups"
print

matrix = group_matrices_one_sample(all_exons, sel_exons, matrix).set_index('name')

#Filter in/out columns in the dataframe

filter_out = False
filter_in = False

if filter_out == False:
    pass
else:
    print "Filtering out columns"
    print
    out_cols = open(str(args.filter_out)).read().split(',')
    out_cols = [w.replace('\n', '') for w in out_cols]
    matrix = matrix.drop(out_cols,1)
    
if filter_in == False:
    pass
else:
    print "Selecting columns"
    print
    in_cols = open(str(args.filter_in)).read().split(',')
    in_cols = [w.replace('\n', '') for w in in_cols]
    matrix = matrix[in_cols]

##Run statistical analysis on the dataset with the "get_statistical_data_for_features" function:

print "Starting statistical analysis"
print

st = get_statistical_data_for_features(matrix, correction_type=args.padj)

print "Finished statistical analysis"
print
##########################################
## Plot CDF for the features

print "Output CDF plots for each features in matrix"
print 

for i in range(len(st)):
    feature = st['Feature'][i]
    name = (st['Feature'][i]).split("/")[-1]
    plt.figure(figsize=(8,8))              
    plot_cdf(matrix[matrix['group'] == 0][str(feature)].values, bins=100, 
             label='Exons in 0', c='black', linewidth=1.5, linestyle='dotted');
    plot_cdf(matrix[matrix['group'] == 1][str(feature)].values, bins=100, 
             label='Exons in 1', c='black', linewidth=1.5, linestyle='dashed');

    plt.legend(loc=0)
    plt.ylim(0,1)

    if feature.find("%") != -1:
        plt.xlim(0,1)
    elif feature.find("phastCon") != -1:
        plt.xlim(0,1)
    elif feature.find("CpG") != -1:
        plt.xlim()
    else:
        plt.xscale('symlog')

    plt.xlabel(str(name), fontsize=14)
    plt.ylabel('Cumulative distribution of samples', fontsize=14)
    plt.title('KS adj. p-val: '+str(round(st['adj_pval_0_vs_1'][i],4)), y=1.01,fontsize=14)
    plt.rc('xtick', labelsize=14)
    plt.rc('ytick', labelsize=14)

    #plt.savefig('/mnt/md0/BioInfo/Project_TARDBP/de_up_exon_gl_con_cum_distribution.pdf', dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
print "Finished CDF plots for features in matrix"
print
