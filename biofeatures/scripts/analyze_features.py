#!/usr/bin/env python
#### Load the required packages

import sys
import argparse
from subprocess import Popen
import pandas as pd
import numpy as np
from pybedtools import BedTool
from scipy import stats
import matplotlib
matplotlib.use('Agg')
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble.partial_dependence import plot_partial_dependence
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
from sklearn.model_selection import cross_val_score
from sklearn.feature_selection import mutual_info_classif,VarianceThreshold
import multiprocessing as mp

pd.options.mode.chained_assignment = None  # default='warn'

plt.ioff()

gandalf = """
    -----------------------       ....
    | YOU SHALL NOT PARSE |     .'' .'''
.   -----------------------   .'   :
\\            \   \         .:    :
 \\             \ \        _:    :       ..----.._
  \\              \     .:::.....:::.. .'         ''.
   \\                 .'  #-. .-######'     #        '.
    \\                 '.##'/ ' ################       :
     \\                  #####################         :
      \\               ..##.-.#### .''''###'.._        :
       \\             :--:########:            '.    .' :
        \\..__...--.. :--:#######.'   '.         '.     :
        :     :  : : '':'-:'':'::        .         '.  .'
        '---'''..: :    ':    '..'''.      '.        :'
           \\  :: : :     '      ''''''.     '.      .:
            \\ ::  : :     '            '.      '      :
             \\::   : :           ....' ..:       '     '.
              \\::  : :    .....####\\ .~~.:.             :
               \\':.:.:.:'#########.===. ~ |.'-.   . '''.. :
                \\    .'  ########## \ \ _.' '. '-.       '''.
                :\\  :     ########   \ \      '.  '-.        :
               :  \\'    '   #### :    \ \      :.    '-.      :
              :  .'\\   :'  :     :     \ \       :      '-.    :
             : .'  .\\  '  :      :     :\ \       :        '.   :
             ::   :  \\'  :.      :     : \ \      :          '. :
             ::. :    \\  : :      :    ;  \ \     :           '.:
              : ':    '\\ :  :     :     :  \:\     :        ..'
                 :    ' \\ :        :     ;  \|      :   .'''
                 '.   '  \\:                         :.''
                  .:..... \\:       :            ..''
                 '._____|'.\\......'''''''.:..'''
                            \\
"""

##Load the parser for arguments

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        print
        self.print_help()
        print
        print(gandalf)
        print    
        print
        print("The following error ocurred in argument parsing:")
        sys.stderr.write('error: %s\n' % message)
        print
        print(
            "Check the help and try to fix the arguments. If the error persists, please contact the corresponding author")
        print   
        sys.exit(2)

class BlankLinesHelpFormatter (argparse.HelpFormatter):
    def _split_lines(self, text, width):
        return super()._split_lines(text, width) + ['']

## Assign input data as system variables

parser = MyParser(description='', 
                  formatter_class=BlankLinesHelpFormatter)

parser.add_argument('-b', '--bed', dest="bed_file",
                    help="BED file with exons/regions of interest.",
                    required=True, metavar='sample.data.bed')

parser.add_argument('-m', '--matrix', dest="matrix",
                    help="Data matrix with biological features created by 'build_datamatrix.py'",
                    required=True, metavar='input.datamatrix.tsv')

parser.add_argument('-o', '--outfile', dest="prefix",  # type=str,
                    help="prefix for use on the output files",
                    metavar="prefix", required=True)

parser.add_argument('-f', '--filter_columns', dest="filter_out",
                    default=False, 
                    help="Text file containing a comma-separated list with names of the columns to be removed from the dataframe in the analysis. Default: False",
                    metavar="filter_out.txt", required=False)

parser.add_argument("-p", '--p_adjust', dest="padj", default='fdr',
                    help="Type of p-value correction used after Kolmogorov-Smirnov test, available options are: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' and 'none'. Default:'fdr'",
                    type=str, metavar='fdr', required=False)

parser.add_argument("-pth", '--p_adjust_threshold', dest="p_th", default=0.05,
                    help="Threshold of adjusted p-value for significance. If using --sig-only-CLF, only significantly different features are passed down to the classifier for group separation. Default: 0.05",
                    type=float, required=False, metavar=0.05)

parser.add_argument('-s','--sig-only', dest="ks_filter", action="store_true",
                    default=False,
                    help="Use only the statistically significant features (found by KS test) in the plotting classification step. Useful for filtering large data matrices to reduce computational time. Can use the '-pth' option to select the threshold of significante for feature selection. Default: False",
                    required=False)

parser.add_argument("-c", '--correlation-filter', dest="corr_th", default=False,
                    help="Remove features which exhibit linear correlation scores above a certain threshold (from 0 to 1). Default: 0.95",
                    type=float, metavar=0.95, required=False)

parser.add_argument("-a", '--ami-filter', dest="ami_th", default=0.01,
                    help="Remove non-informative features based on adjusted mutual information (aMI) scores below a certain threshold (from 0 to 1). If set to 0, will remove features with 0 aMI scores and keep all features with positive scores. Default: 0.01",
                    type=float, metavar=0.01, required=False)

parser.add_argument("-r", '--runs', dest="runs", default=10,
                    help="Number of times (repetitions) to run the classification step. Default:10",
                    type=int, metavar=10)

parser.add_argument("-n", '--sample_size', dest="nsample",
                    type=int, metavar=1, default=1,
                    help="Relative size of randomly sampled exons in comparisson to input exons. Default:1 (i.e. 1x the amount of input exons)")

parser.add_argument("-t", '--train_size', dest="train_size",
                    default=0.80, metavar=0.80,
                    help="Fraction of sample used for training the classifier model. The remaining sample pool will be used for testing the classifier (from 0 to 1). Default: 0.80",
                    type=float)

parser.add_argument("-pr", '--parameters', dest="clf_params",
                    default='optimize',
                    help="Type of parameter selection to be used by the classifier. Available options are: 'optimize' (runs an optimization step with GridSearchCV), 'default' (uses default parameters from GradientBoostClassifier) and 'file' (take in input txt file with a dictionary-like struture with classifier parameters, requires the use of -pf option). Options: 'default', 'optimize' and 'file'. Default:'optimize'",
                    type=str, metavar='optimize', required=False)

parser.add_argument("-scr", '--scoring_metric', dest="scr_metric",
                    default='roc_auc',
                    help="Type of scoring metric used to optimize the classifier hyperparameters if the 'optmize' param option is selected. Possible options include: 'roc_auc', 'accuracy', 'balanced_accuracy', 'precision', 'average_precision' and 'recall' (for a full list of options visit: http://scikit-learn.org/stable/modules/model_evaluation.html). Default: 'roc_auc'",
                    type=str, metavar='roc_auc', required=False)

parser.add_argument("-pf", '--param_file', dest="param_file",
                    help="Input text with with dictionary-like structure with parameter options for GradientBoostClassifier. Ex. {'n_estimators':300,'loss':'deviance',...}",
                    metavar='param_file.txt', required=False)

parser.add_argument("--no-plot", dest="dont_plot_cdf",
                    action="store_true", default=False,
                    help="Use this flag if you want to skip plotting graphs for each feature in the matrix. Default: False")

parser.add_argument("--no-CLF", dest="dont_run_clf",
                    action="store_true", default=False,
                    help="Use this flag if you want to skip the classifying with GradientBoost. Default: False")

parser.add_argument("--ncores", dest="ncores", default=(mp.cpu_count() - 1),
                    help="Number of CPU cores used to multiple jobs on the classifier. Default:(ALL_CORES)-1",
                    type=int, metavar=4)

parser.add_argument('-u','--unstranded', dest="unstranded",
                    action="store_true", default=False, required=False,
                    help="Use this flag if your input file does not contain strand information for genomic intervals. Default: False")

args = parser.parse_args()

if args.unstranded is True:
    strd = False
    print("")
    print("Running in unstranded mode")
    print("")
else:
    strd = True
    print("")
    print("Running in stranded mode")
    print("")

##Define the functions which will be used during the analysis

def group_matrices_one_sample(bt, bt_a, matrix):
    # TODO: refactor common intersect call
    feature_a = bt[0]
    feature_b = bt_a[0]

    if not ((feature_b.strand == "+")  or (feature_b.strand == "-" )) and not ((feature_a.strand == "+")  or (feature_a.strand == "-" )):
        pass
    elif ((feature_b.strand == "+")  or (feature_b.strand == "-" )) and ((feature_a.strand == "+" ) or (feature_a.strand == "-") ):
        pass
    else:
        print("Strand information on input does not match data on matrix. Check your input data.")
        print()
        print("Bed strand data: "+str(feature_b.strand))
        print("Matrix strand data: "+str(feature_a.strand))
        print()
        print("Exiting now. Thanks for using biofeatures!")
        sys.exit()
        
    int_a = bt.intersect(bt_a,
                         s=strd,
                         sorted=True).to_dataframe().drop_duplicates()

    matrix['group'] = 0
    matrix['group'][matrix['name'].isin(int_a['name'])] = 1

    return matrix

def cdf(data, bins=50):
    # data = np.ma.masked_array(data, np.isnan(data))
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
    ax.plot(x, y, **plotting_args)

def get_mean_and_std(df):
    df2 = df.T
    cols = df2.columns
    df2['mean'] = df2.apply(lambda x: np.mean(x[cols]), 1)
    df2['std'] = df2.apply(lambda x: np.std(x[cols]), 1)
    return df2

def plot_barcharts(df, title, save):
    df_t = df.T
    cols = df_t.columns

    df_t['mean'] = df_t.apply(lambda x: np.mean(x[cols]), 1)
    df_t['std'] = df_t.apply(lambda x: np.std(x[cols]), 1)

    df_t.to_csv('./' + args.prefix + '.analysis/classifier_metrics/' + save + '.tsv', sep='\t')

    N = df_t.shape[0]
    means = df_t['mean']
    std = df_t['std']

    ind = np.arange(N)  # the x locations for the groups
    width = 0.5  # the width of the bars

    fig, ax = plt.subplots()
    ax.bar(ind, means, width, color='lightgrey', yerr=std,
                    ecolor='black')

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Value', fontsize=14)
    ax.set_title(title, fontsize=14)
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(df_t.index, rotation='vertical')


def plot_barchart_importance(df):
    importance_ind = df.set_index("Feature")

    N = 10
    ind = np.arange(N)  # the x locations for the groups
    width = 0.5

    rel = \
        importance_ind.sort_values('mean_rel_importance',
                                   ascending=False)[
            'mean_rel_importance'].head(N)
    rel_err = \
        importance_ind.sort_values('mean_rel_importance',
                                   ascending=False)[
            'std_rel_importance'].head(N)

    fig, ax = plt.subplots()
    ax.bar(ind, rel, width, color='lightgrey', yerr=rel_err, ecolor='black')
    ax.set_ylabel('Relative importance', fontsize=14)
    ax.set_title('Mean relative importance values', fontsize=14)
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(rel.index, rotation='vertical', fontsize=14)
    plt.savefig(
        './' + args.prefix + '.analysis/classifier_metrics/mean_relative_importance.pdf',
        dpi=300, bbox_inches='tight')
    plt.close()

    raw = \
        importance_ind.sort_values(by=['mean_raw_importance'],
                                   ascending=False)[
            'mean_raw_importance'].head(N)
    raw_err = \
        importance_ind.sort_values(by=['mean_raw_importance'],
                                   ascending=False)[
            'std_raw_importance'].head(N)

    fig, ax = plt.subplots()
    ax.bar(ind, raw, width, color='lightgrey', yerr=raw_err, ecolor='black')
    ax.set_ylabel('Importance', fontsize=14)
    ax.set_title('Mean importance values', fontsize=14)
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(raw.index, rotation='vertical', fontsize=14)
    plt.savefig(
        './' + args.prefix + '.analysis/classifier_metrics/mean_importance.pdf',
        dpi=300, bbox_inches='tight')
    plt.close()


##Create directory for output

Popen('mkdir -p ./' + args.prefix + '.analysis', shell=True)
Popen('mkdir -p ./' + args.prefix + '.analysis/feature_plots', shell=True)
Popen('mkdir -p ./' + args.prefix + '.analysis/classifier_metrics', shell=True)
Popen('mkdir -p ./' + args.prefix + '.analysis/statistical_analysis', shell=True)

    
##Load the bed file created with all exons

print("Loading bed file with regions of interest")
print()


bed_input = BedTool(args.bed_file).sort()
    
##Load the datamatrix generatade by "buildadatamatrix.py"

print("Loading datamatrix")
print()

matrix = pd.concat(pd.read_table(args.matrix, iterator=True, chunksize=10000),
                   ignore_index=True).set_index('name').drop_duplicates().reset_index()

# Filter in/out columns in the dataframe

if not args.filter_out:
    pass
else:
    print("Filtering out columns")
    print()
    out_cols = open(str(args.filter_out)).read().split(',')
    out_cols = [w.replace('\n', '') for w in out_cols]
    matrix = matrix.drop(out_cols, 1)

##Intersect the exons found in the analysis to get groups 1 (positive) and 0 (negative) in the matrix

print("Finding input exons in the matrix and selecting groups")
print()

bed_from_matrix = pd.DataFrame()
bed_from_matrix['chr'] = matrix['name'].apply(lambda x: x.split('_')[3],1)
bed_from_matrix['start'] = matrix['name'].apply(lambda x: x.split('_')[4],1)
bed_from_matrix['end'] = matrix['name'].apply(lambda x: x.split('_')[5],1)
bed_from_matrix['name'] = matrix['name']
bed_from_matrix['score'] = 0

if strd is True:
    bed_from_matrix['strand'] = matrix['name'].apply(lambda x: x.split('_')[6],1)
else:
    pass

bed_from_matrix = BedTool.from_dataframe(bed_from_matrix).sort()

matrix = group_matrices_one_sample(bed_from_matrix, bed_input, matrix).set_index('name')
   
print("Starting statistical analysis")
print()
print("Calculating Komlogorov-Smirnov test for each feature in the matrix")
print()

features = list(matrix.drop('group',1).columns)
df_list = []

for i in range(len(features)):
    sl = matrix[[features[i],'group']]
    res = stats.ks_2samp(sl[sl['group'] == 0][features[i]].astype(float),
                         sl[sl['group'] == 1][features[i]].astype(float))
    d = {'Feature': features[i], 'ks': res[0], 'pval': res[1]}
    df = pd.DataFrame(data=d, index=np.arange(1))
    df_list.append(df)

st = pd.concat(df_list, 0).reset_index().drop('index',1)
    
print("Adjusting pvalues using "+str(args.padj)+" and saving output")
print()
    
statsR = importr('stats')
st['adj_pval'] = statsR.p_adjust(FloatVector(st['pval']),method=str(args.padj))
st.to_csv('./' + args.prefix + '.analysis/statistical_analysis_output.tsv',
          sep='\t', index=False)
st.to_excel('./' + args.prefix + '.analysis/statistical_analysis_output.xlsx',
            index=False)

print("Finished statistical analysis")
print()
    
if not args.ks_filter:
    pass
elif args.ks_filter:
    print("Filtering statistically significantly features.")
    print()
    sig_only = st[st['adj_pval'] <= float(args.p_th)]['Feature'].tolist()
    sig_only.append('group')
    matrix = matrix[sig_only]
    
if not args.corr_th:
    pass
elif args.corr_th:
    print("Calculating linear correlation matrix and filtering features above threshold.")
    print()
    corr_matrix = matrix.drop(['group'], axis=1).corr().abs()
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(np.bool))
    to_drop = [column for column in upper.columns if any(upper[column] >= args.corr_th)]
    matrix = matrix.drop(to_drop, axis=1)
    
    print('Features above threshold:')
    print()
    print(to_drop)
    print()
    print('Saving correlation matrix')
    print()
    corr_matrix.to_csv('./' + args.prefix + '.analysis/statistical_analysis/correlation_matrix.tsv',
                       sep='\t')
    corr_matrix.to_excel('./' + args.prefix + '.analysis/statistical_analysis/correlation_matrix.xlsx')
    
    
if not args.dont_plot_cdf:
    print("Output plots for each features in matrix")
    print()
    features = list(matrix.drop('group',1).columns)
    
    title_size=16
    tick_size=14
    axis_size=14
    
    for i in range(len(features)):
        name = (features[i]).split("/")[-1]
        sl = matrix[[features[i],'group']]

        plt.figure(1, figsize=(18,6))
        plt.suptitle("Cumulative distribution and kernel density for "+str(name), 
                     fontsize=title_size, y=1.01)
        plt.subplot(121)

        sns.set(font_scale = 1.5)
        sns.set_style("whitegrid", {'axes.grid' : False})
    #     plt.figure(figsize=(8, 8))
        plot_cdf(sl[sl['group'] == 0][features[i]].values,
                 bins=100,
                 label='Background regions', c='black', linewidth=1.5,
                 linestyle='solid')
        plot_cdf(sl[sl['group'] == 1][features[i]].values,
                 bins=100,
                 label='Input regions', c='black', linewidth=1.5,
                 linestyle='dashed')

        plt.legend(loc=0)
        plt.ylim(0, 1)

        if name.find("%") != -1:
            plt.xlim(0, 1)
        elif name.find("phastCon") != -1:
            plt.xlim(0, 1)
#        elif name.find("_count_") != -1:
#            plt.xscale('symlog')
        else:
            plt.xlim()

        plt.xlabel(str(name), fontsize=axis_size)
        plt.ylabel('Cumulative distribution of samples', fontsize=axis_size)
        plt.rc('xtick', labelsize=tick_size)
        plt.rc('ytick', labelsize=tick_size)

        plt.subplot(122)
        #sns.set(font_scale = 1.5)
        sns.set_style("whitegrid", {'axes.grid' : False})

        shade = False
        label0='Background regions'
        label1='Input regions'

        if name.find('_count_') != -1 or \
           name.find('QGRS') != -1 or \
            name.find('MFE') != -1:
                sns.kdeplot(sl[sl['group'] == 0][features[i]], 
                            color='k', 
                            cut=0, 
                            shade=shade,
                            bw=1,
                            label=label0)

                sns.kdeplot(sl[sl['group'] == 1][features[i]], 
                            color='k', 
                            cut=0, 
                            shade=shade,
                            bw=1,
                            label=label1,
                            ls='--')
        else:
            sns.kdeplot(sl[sl['group'] == 0][features[i]], 
                        color='k', 
                        cut=0, 
                        shade=shade,
                        label=label0)

            sns.kdeplot(sl[sl['group'] == 1][features[i]], 
                        color='k', 
                        cut=0, 
                        shade=shade,
                        label=label1,
                        ls='--')


        plt.xlabel(str(name), fontsize=axis_size)
        plt.ylabel('Density', fontsize=axis_size)
        plt.rc('xtick', labelsize=tick_size)
        plt.rc('ytick', labelsize=tick_size)

        if name.find("%") != -1:
            plt.xlim(0, 1)
        elif name.find("phastCon") != -1:
            plt.xlim(0, 1)
#        elif name.find("_count_") != -1:
#            plt.xscale('symlog')
        else:
            plt.xlim()
        plt.savefig('./' + args.prefix + '.analysis/feature_plots/' + name + '.pdf',
                    dpi=300, bbox_inches='tight')
        plt.close()

    print("Finished CDF plots for features in matrix") 
    print()
elif args.dont_plot_cdf:
    pass

if args.dont_run_clf:
    print("Analysis complete. Thank you for using biofeatures.")
    print()
    sys.exit()
else:
    pass

##For the classification step, there are a few issues that we need to solve:
# 1. In order to get an accurate background representation, we need to sample a large portion of the exons
# 2. On the other hand, using a higher ammount of exons can introduce a lot of false positives. This leads 
# to a loss in the positive predictive power (or Precision) of the classifier.
# 3. So we need to make a trade-off for the accuracy of the biological features and the precision of classification
# 4. Our approach is to use a small subset of data for background in classification (3x the size of the test set) 
# and perform multiple classifications with random small backgrouds. The feature importance output is the median
# of the importance in each run with the associated standard deviation.

runs = np.arange(args.runs)
sample_size = args.nsample
gbclf_params = str(args.clf_params)
train_size = args.train_size
scr_metric = args.scr_metric
ncores = args.ncores

##Popen('mkdir -p ./' + args.prefix + '.analysis/classifier_metrics', shell=True)

importance = st.copy()
    
deviance_train = pd.DataFrame()
deviance_test = pd.DataFrame()

msr = pd.DataFrame(columns=['aMI', 'MSE'])

mtc = pd.DataFrame(
    columns=['Accuracy', 'P.P.V.', 'N.P.V.', 'Sensitivity', 'Specificity'])

conf_df = pd.DataFrame(columns=['T.N.', 'T.P.', 'F.N.', 'F.P.'])

roc_tpr_all = pd.DataFrame()
roc_fpr_all = pd.DataFrame()
roc_auc_all = pd.DataFrame()

pre_all = pd.DataFrame()
rec_all = pd.DataFrame()
av_pre_all = pd.DataFrame()

if gbclf_params == 'default':
    print("Using default parameters.")
    print()

if gbclf_params == 'file':
    print("Using input file as param dictionary")
    print()
    bp = eval(open(str(args.param_file)).read())
    
if gbclf_params == 'optimize':
    
    print("Starting parameter tuning process with stratified 10-fold cross validation")
    print()
    print("Selected scoring parameter: "+str(scr_metric))
    print()

    ##Since we are using the GradientBooost Classfier, the first step is to find how many features are 
    ##optimal for group classification. To do that, we use recursive feature elimination with cross-validation.
    ##For this step, we will use the StratifiedKFold approach with 10-fold cross validation. The "accuracy" 
    ##scoring is proportional to the number of correct classifications

    
    df_cl = matrix[matrix['group'] != 0]
    df_zero = matrix[matrix['group'] == 0].sample(df_cl.shape[0] * sample_size)
    Z = pd.concat([df_cl, df_zero]).drop_duplicates()
    
    y = np.array(Z.group).astype('int')
    X = np.array(Z.drop('group', 1))
    
    print('Removing features non-informative features based on aMI score')
    print()
    
    mi_scores = mutual_info_classif(X, y, 
                            discrete_features='auto', 
                            n_neighbors=3, 
                            copy=True, 
                            random_state=None)

    dd = []

    for score, fname in sorted(zip(mi_scores, Z.columns.tolist()), reverse=True):
        dd.append([fname, score])

    dm = pd.DataFrame(dd, columns=['features','ami_score'])

    selector = VarianceThreshold()
    selector.fit(X, y)
    var_scores = list(selector.variances_)

    dd = []

    for score, fname in sorted(zip(var_scores, Z.drop('group', axis=1).columns.tolist()), reverse=True):
        dd.append([fname, score])

    dv = pd.DataFrame(dd, columns=['features','var_score'])

    ami_var = dm.merge(dv, on='features')
    
    ami_var.to_csv('./' + args.prefix + '.analysis/statistical_analysis/features_ami_and_variance_scores.tsv',
                       sep='\t')
    ami_var.to_excel('./' + args.prefix + '.analysis/statistical_analysis/features_ami_and_variance_scores.xlsx')
    
    if args.ami_th == 0:
        if len(ami_var[ami_var['ami_score'] == args.ami_th]['features'].tolist()):
            Z = Z.drop(ami_var[ami_var['ami_score'] == args.ami_th]['features'].tolist(), axis=1)
            y = np.array(Z.group).astype('int')
            X = np.array(Z.drop('group', 1))
        else:
            print('No features found with aMI lower than threshold')
            print()
            pass
    else:
        if len(ami_var[ami_var['ami_score'] < args.ami_th]['features'].tolist()):
            Z = Z.drop(ami_var[ami_var['ami_score'] < args.ami_th]['features'].tolist(), axis=1)
            y = np.array(Z.group).astype('int')
            X = np.array(Z.drop('group', 1))
        else:
            print('No features found with aMI lower than threshold')
            print()
            pass
        
    X_train, X_test, y_train, y_test = train_test_split(X,
                                                        y,
                                                        train_size=train_size,
                                                        test_size=(1-train_size))
    
    print('Optmizing number of boosting stages (estimators)')
    print()

    opt = GradientBoostingClassifier(warm_start=True,
                                     learning_rate=0.1, 
                                     loss='deviance',
                                     subsample=0.8,
                                     random_state=0,
                                     max_depth=6, 
                                     min_samples_split=0.01, 
                                     min_samples_leaf=0.001,
                                     max_features='sqrt',
                                     validation_fraction=0.2,
                                     n_iter_no_change=5, 
                                     tol=1e-4)

    grid = {'n_estimators': range(10,301,10)}

    opt_s1 = GridSearchCV(estimator=opt, param_grid=grid, n_jobs=ncores,
                          cv=StratifiedKFold(10), scoring=scr_metric, return_train_score=True)

    opt_s1.fit(X, y)

    print("SCORES:")
    print(pd.DataFrame(opt_s1.cv_results_)[['params','mean_test_score','std_test_score']])
    print()
    print('Best result:')
    print('SCORE: '+str(opt_s1.best_score_)+'; PARAM: '+str(opt_s1.best_params_)+'; ESTIMATORS USED WITH EARLY STOPPING: '+str(opt_s1.best_estimator_.n_estimators_))
    print()

    n_est = list(opt_s1.best_params_.values())[0]
    
    print('Optmizing number of maximum depth of the individual estimators and minimum number of samples to split a node')
    print()

    opt = GradientBoostingClassifier(warm_start=True,
                                     n_estimators=n_est,
                                     learning_rate=0.1, 
                                     loss='deviance',
                                     subsample=0.8,
                                     random_state=0,
                                     min_samples_leaf=0.001,
                                     max_features='sqrt',
                                     validation_fraction=0.2,
                                     n_iter_no_change=5, 
                                     tol=1e-4)

    grid = {'max_depth':range(2,21,1), 
            'min_samples_split':list(np.arange(0.001, 0.021, 0.001).round(3))}

    opt_s2 = GridSearchCV(estimator=opt, param_grid=grid, n_jobs=ncores,
                          cv=StratifiedKFold(10), scoring=scr_metric, return_train_score=True)

    opt_s2.fit(X, y)

    print("SCORES:")
    print(pd.DataFrame(opt_s2.cv_results_)[['params','mean_test_score','std_test_score']])
    print()
    print('Best result:')
    print('SCORE: '+str(opt_s2.best_score_)+'; PARAM: '+str(opt_s2.best_params_))
    print()
    
    m_depth = list(opt_s2.best_params_.values())[0]
    n_split = list(opt_s2.best_params_.values())[1]
    
    print('Adjusting number of minimum samples in each leaf after split')
    print()

    opt = GradientBoostingClassifier(warm_start=True,
                                     n_estimators=n_est,
                                     min_samples_split=n_split,
                                     max_depth=m_depth,
                                     learning_rate=0.1, 
                                     loss='deviance',
                                     subsample=0.8,
                                     random_state=0,
                                     max_features='sqrt',
                                     validation_fraction=0.2,
                                     n_iter_no_change=5, 
                                     tol=1e-4)

    grid = {'min_samples_split':list(np.array(np.linspace(n_split/2, n_split*2, 10).round(5))),
            'min_samples_leaf':list(np.array(np.linspace(n_split/20, n_split*2, 10).round(5)))}

    opt_s3 = GridSearchCV(estimator=opt, param_grid=grid, n_jobs=ncores,
                          cv=StratifiedKFold(10), scoring=scr_metric, return_train_score=True)

    opt_s3.fit(X, y)

    print("SCORES:")
    print(pd.DataFrame(opt_s3.cv_results_)[['params','mean_test_score','std_test_score']])
    print()
    print('Best result:')
    print('SCORE: '+str(opt_s3.best_score_)+'; PARAM: '+str(opt_s3.best_params_))
    print()
    
    n_leaf = list(opt_s3.best_params_.values())[0]
    n_split = list(opt_s3.best_params_.values())[1]

    print('Optmizing maximmum number of features')
    print()

    sqrt = int(np.sqrt(X.shape[1]))

    opt = GradientBoostingClassifier(warm_start=True,
                                     n_estimators=n_est,
                                     min_samples_split=n_split,
                                     min_samples_leaf=n_leaf,
                                     max_depth=m_depth,
                                     learning_rate=0.1, 
                                     loss='deviance',
                                     subsample=0.8,
                                     random_state=0,
                                     validation_fraction=0.2,
                                     n_iter_no_change=5, 
                                     tol=1e-4)                                 

    grid = {'max_features':list(np.arange(((sqrt/2)), ((sqrt*2)+2), 2, dtype=int))}

    opt_s4 = GridSearchCV(estimator=opt, param_grid=grid, n_jobs=ncores,
                          cv=StratifiedKFold(10), scoring=scr_metric, return_train_score=True)

    opt_s4.fit(X, y)

    print("SCORES:")
    print(pd.DataFrame(opt_s4.cv_results_)[['params','mean_test_score','std_test_score']])
    print()
    print('Best result:')
    print('SCORE: '+str(opt_s4.best_score_)+'; PARAM: '+str(opt_s4.best_params_))
    print()
    
    n_features = list(opt_s4.best_params_.values())[0]
    
    print('Adjusting number of estimators and learning rate')
    print()

    opt = GradientBoostingClassifier(warm_start=True,
                                     min_samples_split=n_split,
                                     min_samples_leaf=n_leaf,
                                     max_depth=m_depth,
                                     loss='deviance',
                                     subsample=0.8,
                                     random_state=0,
                                     max_features=n_features,
                                     validation_fraction=0.2,
                                     n_iter_no_change=5, 
                                     tol=1e-4)

    grid = {'n_estimators':[n_est, n_est*5, n_est*10, n_est*20],
            'learning_rate':[0.1, 0.02, 0.01, 0.005]}

    opt_s5 = GridSearchCV(estimator=opt, param_grid=grid, n_jobs=ncores,
                          cv=StratifiedKFold(10), scoring=scr_metric, return_train_score=True)

    opt_s5.fit(X, y)

    print("SCORES:")
    print(pd.DataFrame(opt_s5.cv_results_)[['params','mean_test_score','std_test_score']])
    print()
    print('Best result:')
    print('SCORE: '+str(opt_s5.best_score_)+'; PARAM: '+str(opt_s5.best_params_))
    print()
    
    n_est = list(opt_s5.best_params_.values())[1]
    l_rate = list(opt_s5.best_params_.values())[0]
    
    print('Scoring the optimized hyperparameters with stratified 10-fold cross-validation')

    scr = GradientBoostingClassifier(warm_start=True,
                                     n_estimators=n_est,
                                     min_samples_split=n_split,
                                     min_samples_leaf=n_leaf,
                                     max_depth=m_depth,
                                     learning_rate=l_rate, 
                                     loss='deviance',
                                     subsample=0.8,
                                     random_state=0,
                                     max_features=n_features,
                                     validation_fraction=0.2,
                                     n_iter_no_change=5, 
                                     tol=1e-4).fit(X_train, y_train)

    scores = cross_val_score(scr, X_test, y_test, cv=StratifiedKFold(10), scoring=scr_metric)
    print()
    print("Final 10-fold CV score for optimized hyperparameters: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
    print()
    
    bp = scr.get_params()
    
    print("Exporting final hyperparameter selection to 'best_params.txt' file:")
    print()
    print(bp)
    
    with open('./' + args.prefix \
      + '.analysis/best_params.txt','w') as fp:
        fp.write(str(scr.get_params()))

for run_i in range(len(runs)):
    run_id = str(runs[run_i] + 1)
    df_cl = matrix[matrix['group'] != 0]
    df_zero = matrix[matrix['group'] == 0].sample(df_cl.shape[0] * sample_size)
    Z = pd.concat([df_cl, df_zero]).drop_duplicates()

    Popen(
        'mkdir -p ./' + args.prefix + '.analysis/classifier_metrics/run_' + run_id,
        shell=True)

    print(("Run ID: " + str(run_id) + ', Input size: ' + str(
        df_cl.shape[0]) + ', Background size: ' + str(df_zero.shape[0])))
    print()
    print(
        "Starting classification analysis with GradientBoost. This may take a long time depending on the size of the dataset")
    print()

    ##The first step is separating the groups from the data in different arrays

    print("Loading data in arrays")
    print()

    y = np.array(Z.group).astype('int')
    X = np.array(Z.drop('group', 1))
    
    print('Removing non-informative features based on aMI score')
    print()
    
    mi_scores = mutual_info_classif(X, y, 
                            discrete_features='auto', 
                            n_neighbors=3, 
                            copy=True, 
                            random_state=None)

    dd = []

    for score, fname in sorted(zip(mi_scores, Z.columns.tolist()), reverse=True):
        dd.append([fname, score])

    dm = pd.DataFrame(dd, columns=['features','ami_score'])

    selector = VarianceThreshold()
    selector.fit(X, y)
    var_scores = list(selector.variances_)

    dd = []

    for score, fname in sorted(zip(var_scores, Z.drop('group', axis=1).columns.tolist()), reverse=True):
        dd.append([fname, score])

    dv = pd.DataFrame(dd, columns=['features','var_score'])

    ami_var = dm.merge(dv, on='features')
    
    ami_var.to_csv('./' + args.prefix + '.analysis/classifier_metrics/run_' + run_id + '/run_' + run_id + '_features_ami_and_variance_scores.tsv',
                       sep='\t')
    ami_var.to_excel('./' + args.prefix + '.analysis/classifier_metrics/run_' + run_id + '/run_' + run_id + '_features_ami_and_variance_scores.xlsx')
    
    if args.ami_th == 0:
        if len(ami_var[ami_var['ami_score'] == args.ami_th]['features'].tolist()) > 0:
            Z = Z.drop(ami_var[ami_var['ami_score'] == args.ami_th]['features'].tolist(), axis=1)
            y = np.array(Z.group).astype('int')
            X = np.array(Z.drop('group', 1))
        else:
            print('No features found with aMI lower than threshold')
            print()
            pass
    else:
        if len(ami_var[ami_var['ami_score'] < args.ami_th]['features'].tolist()) > 0:
            Z = Z.drop(ami_var[ami_var['ami_score'] < args.ami_th]['features'].tolist(), axis=1)
            y = np.array(Z.group).astype('int')
            X = np.array(Z.drop('group', 1))
        else:
            print('No features found with aMI lower than threshold')
            print()
            pass
        
    ##Then, we split the total data in 2 groups: training group (which we will train our model and optimize the parameters) 
    ## and the test group, which we will use for evaluating the performance of the classifier. The default test size is 0.25.

    print("Splitting train and test datasets")
    print()

    X_train, X_test, y_train, y_test = train_test_split(X,
                                                        y,
                                                        train_size=train_size,
                                                        test_size=(1-train_size))

    
    print("Fitting the GradientBoost model with the training set")
    print()

    names = Z.drop('group', 1).columns

    if gbclf_params == 'default':
        clf = GradientBoostingClassifier(warm_start=True,
                                         validation_fraction=0.2,
                                         n_iter_no_change=5, 
                                         tol=1e-4)
    else:
        clf = GradientBoostingClassifier(**bp)

    clf.fit(X_train, y_train)

    print((
        "Extracting feature importance for run " + run_id + " and merging with statistical data"))
    print()

    clf_fi = []

    for score, fname in sorted(zip(clf.feature_importances_, 
                                   Z.drop('group', axis=1).columns.tolist()), 
                               reverse=True):
            clf_fi.append([fname, score])
        
    a = pd.DataFrame(clf_fi, 
                    columns=['Feature','run_' + run_id + '_raw_importance']
                    ).sort_values('run_' + run_id + '_raw_importance',
                                  ascending=False)

    max_importance = a['run_' + run_id + '_raw_importance'].max()
    a['run_' + run_id + '_rel_importance'] = 100.0 * (
        a['run_' + run_id + '_raw_importance'] / max_importance)
    
    importance = importance.merge(a, on='Feature')

    print("Plotting partial dependance and deviance scores from classifier")
    print()

    # Select top 5 most important features and plot partial dependance and interaction plots

    do_these = a.sort_values(by=('run_' + run_id + '_raw_importance'), 
                             ascending=False).index[:10].values

    names = a.sort_values(by=('run_' + run_id + '_raw_importance'), 
                      ascending=False)[:10]['Feature']

    features = list(do_these)
    features.append((do_these[0], do_these[1]))
    fig, axs = plot_partial_dependence(clf, X_train, features,
                                       feature_names=names,
                                       n_jobs=ncores, grid_resolution=100,
                                       figsize=(12,16),
                                       label=0)

    plt.subplots_adjust(top=0.9)  # tight_layout causes overlap with suptitle

    plt.tight_layout()
    plt.savefig(
        './' + args.prefix + '.analysis/classifier_metrics/run_' + run_id + '/run_' + run_id + '_partial_dependance.pdf',
        dpi=600,
        bbox_inches='tight')
    plt.close()

    ##Plot deviance x number os iteration
    ###############################################################################
    test_score = np.zeros((clf.estimators_.shape[0],), dtype=np.float64)

    for i, y_pred in enumerate(clf.staged_decision_function(X_test)):
        test_score[i] = clf.loss_(y_test, y_pred)

    plt.figure(figsize=(6, 6))
    plt.subplot(1, 1, 1)
    plt.title('Deviance')
    plt.plot(np.arange(clf.estimators_.shape[0]) + 1, clf.train_score_, 'b-',
             label='Training Set Deviance')
    plt.plot(np.arange(clf.estimators_.shape[0]) + 1, test_score, 'r-',
             label='Test Set Deviance')
    plt.legend(loc='upper right')
    plt.xlabel('Boosting Iterations')
    plt.ylabel('Deviance')
    plt.savefig(
        '' + args.prefix + '.analysis/classifier_metrics/run_' + run_id + '/run_' + run_id + '_deviance.pdf',
        dpi=600,
        bbox_inches='tight')
    plt.close()

    deviance_train['run_' + run_id] = pd.Series(clf.train_score_)
    deviance_test['run_' + run_id] = pd.Series(test_score)

    print("Run trained classifier in test subset from data")
    print
    predicted_values = clf.predict(X_test)
    mse = mean_squared_error(y_test, predicted_values)
    ami = adjusted_mutual_info_score(y_test, predicted_values)
    conf = confusion_matrix(y_test, predicted_values)
    spe = round((float(conf[0, 0]) / (conf[1, 0] + conf[0, 0])), 2)
    sen = round((float(conf[1, 1]) / (conf[0, 1] + conf[1, 1])), 2)
    npv = round((float(conf[0, 0]) / (conf[0, 0] + conf[0, 1])), 2)
    ppv = round((float(conf[1, 1]) / (conf[1, 0] + conf[1, 1])), 2)
    acc = round((float(conf[0, 0] + conf[1, 1]) / (
        conf[0, 0] + conf[1, 0] + conf[0, 1] + conf[1, 1])), 2)

    msr_r = pd.DataFrame(data=[[ami, mse]], columns=['aMI', 'MSE'])

    msr = msr.append(msr_r).reset_index().drop('index', 1)

    mtc_r = pd.DataFrame(data=[[acc, ppv, npv, sen, spe]],
                         columns=['Accuracy', 'P.P.V.', 'N.P.V.',
                                  'Sensitivity', 'Specificity'])

    mtc = mtc.append(mtc_r).reset_index().drop('index', 1)

    conf_df_r = pd.DataFrame(
        data=[[conf[0, 0], conf[1, 1], conf[1, 0], conf[0, 1]]],
        columns=['T.N.', 'T.P.', 'F.N.', 'F.P.'])

    conf_df = conf_df.append(conf_df_r).astype(int).reset_index().drop('index',
                                                                       1)

    ##For creating ROC and precision curves, we need to binarize the output and get score functions
    print("Binarizing output")
    print()
    y_class = label_binarize(y_test, classes=[0, 1])
    n_classes = y_class.shape[1]
    y_score = clf.fit(X_train, y_train).decision_function(X_test)

    ###Now, let's calculate the ROC score for the classification

    print("Computing ROC score for classifier")
    print()

    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    # Compute micro-average ROC curve and ROC area
    fpr["micro_run_" + run_id], tpr["micro_run_" + run_id], _ = roc_curve(
        y_test.ravel(), y_score.ravel())
    roc_auc["micro_run_" + run_id] = auc(fpr["micro_run_" + run_id],
                                         tpr["micro_run_" + run_id])

    fpr_r = pd.DataFrame(data=fpr["micro_run_" + run_id], columns=['F.P.R.'])
    fpr_r.to_csv(
        './' + args.prefix + '.analysis/classifier_metrics/run_' + run_id + '/run_' + run_id + '_roc_fpr.tsv',
        sep='\t')

    tpr_r = pd.DataFrame(data=tpr["micro_run_" + run_id], columns=['T.P.R.'])
    tpr_r.to_csv(
        './' + args.prefix + '.analysis/classifier_metrics/run_' + run_id + '/run_' + run_id + '_roc_tpr.tsv',
        sep='\t')

    roc_fpr_all = pd.concat([roc_fpr_all, fpr_r], ignore_index=True,
                            axis=1)  # .fillna(1)
    roc_tpr_all = pd.concat([roc_tpr_all, tpr_r], ignore_index=True,
                            axis=1)  # .fillna(1)
    roc_auc_all['run_' + run_id] = np.array([roc_auc['micro_run_' + run_id]])

    # Plot ROC curve

    plt.plot(fpr["micro_run_" + run_id], tpr["micro_run_" + run_id],
             label='ROC curve (AUC = {0:0.2f})'
                   ''.format(roc_auc["micro_run_" + run_id]), linewidth=2)

    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate', size=20)
    plt.ylabel('True Positive Rate', size=20)
    plt.title('Receiver operating characteristic', size=20)
    plt.legend(loc="lower right", fontsize=18)
    plt.savefig(
        './' + args.prefix + '.analysis/classifier_metrics/run_' + run_id + '/run_' + run_id + '_roc_curve.pdf',
        dpi=300,
        bbox_inches='tight')
    plt.close()

    print("Plotting Precision-Recall curve")
    print()
    # Compute Precision-Recall and plot curve
    precision = dict()
    recall = dict()
    average_precision = dict()

    precision["micro_run_" + run_id], recall[
        "micro_run_" + run_id], _ = precision_recall_curve(y_test, y_score)
    average_precision["micro_run_" + run_id] = average_precision_score(y_test,
                                                                       y_score)

    pre_r = pd.DataFrame(data=precision["micro_run_" + run_id],
                         columns=['Precision'])
    pre_r.to_csv(
        './' + args.prefix + '.analysis/classifier_metrics/run_' + run_id + '/run_' + run_id + '_precision.tsv',
        sep='\t')

    rec_r = pd.DataFrame(data=recall["micro_run_" + run_id],
                         columns=['Recall'])
    rec_r.to_csv(
        './' + args.prefix + '.analysis/classifier_metrics/run_' + run_id + '/run_' + run_id + '_recall.tsv',
        sep='\t')

    pre_all = pd.concat([pre_all, pre_r], ignore_index=True,
                        axis=1)  # .fillna(0)
    rec_all = pd.concat([rec_all, rec_r], ignore_index=True,
                        axis=1)  # .fillna(0)
    av_pre_all['run_' + run_id] = np.array(
        [average_precision['micro_run_' + run_id]])

    # Plot Precision-Recall curve

    plt.figure(figsize=(6, 6))

    plt.clf()
    plt.plot(recall["micro_run_" + run_id], precision["micro_run_" + run_id],
             lw=2, color='navy',
             label='Precision-Recall curve (AUC={0:0.2f})'.format(
                 average_precision["micro_run_" + run_id]))
    plt.title('Precision-Recall', size=20)
    plt.xlabel('Recall', size=20)
    plt.ylabel('Precision', size=20)
    plt.xlim([-0.05, 1.05])
    # plt.ylim([-0.05, 1.05])
    plt.ylim(ymax=1.05)
    plt.legend(loc="lower left", fontsize=18)
    plt.savefig(
        './' + args.prefix + '.analysis/classifier_metrics/run_' + run_id + '/run_' + run_id + '_precision_recall_curve.pdf',
        dpi=300,
        bbox_inches='tight')
    plt.close()

    print('=======================')
    print()
    print("Overall classifier scores:")
    print(("Mean Squared Error (MSE): %.4f" % mse))
    print(("adjusted Mutual Information (aMI): %.4f" % ami))
    print(('ROC AUC: {0:0.4f}'.format(roc_auc["micro_run_" + run_id])))
    print(('Precision-Recall AUC: {0:0.4f}'.format(average_precision["micro_run_" + run_id])))
    print()
    print("Confusion matrix")
    print(conf)
    print()
    print(('Accuracy: ' + str(acc*100) + '%'))
    print(('Positive predictive value (Precision): ' + str(ppv*100) + '%'))
    print(('Negative predictive value: ' + str(npv*100) + '%'))
    print(('Sensitivity (Recall): ' + str(sen*100) + '%'))
    print(('Specificity: ' + str(spe*100) + '%'))
    print()

print("Done with classification steps")
print()
print("Calculating median importance values and standard deviation")
print()

cols = importance.filter(like='_raw_importance', axis=1).columns
importance['mean_raw_importance'] = importance.apply(
    lambda x: np.mean(x[cols]), 1)
importance['std_raw_importance'] = importance.apply(lambda x: np.std(x[cols]),
                                                    1)

cols = importance.filter(like='_rel_importance', axis=1).columns
importance['mean_rel_importance'] = importance.apply(
    lambda x: np.mean(x[cols]), 1)
importance['std_rel_importance'] = importance.apply(lambda x: np.std(x[cols]),
                                                    1)
importance = importance[importance.columns.drop(list(importance.filter(regex='run_')))]

importance.to_csv('./' + args.prefix + '.analysis/classifier_importance.tsv',
                  sep='\t')

importance.to_excel('./' + args.prefix + '.analysis/classifier_importance.xlsx')

print("Plotting raw importance and relative importance barcharts")
print()

plot_barchart_importance(importance)

print("Plotting mean deviance curve and error region")
print()
cols = deviance_test.columns
deviance_test['mean'] = deviance_test.apply(lambda x: np.mean(x[cols]), 1)
deviance_test['std'] = deviance_test.apply(lambda x: np.std(x[cols]), 1)

deviance_test.to_csv(
    './' + args.prefix + '.analysis/classifier_metrics/classifier_test_deviance.tsv', sep='\t')

cols = deviance_train.columns
deviance_train['mean'] = deviance_train.apply(lambda x: np.mean(x[cols]), 1)
deviance_train['std'] = deviance_train.apply(lambda x: np.std(x[cols]), 1)

deviance_train.to_csv(
    './' + args.prefix + '.analysis/classifier_metrics/classifier_train_deviance.tsv', sep='\t')

plt.figure(figsize=(6, 6))
plt.subplot(1, 1, 1)
plt.title('Deviance')

train_err = np.array(deviance_train.dropna()['std'])
test_err = np.array(deviance_test.dropna()['std'])

plt.plot(np.arange(len(deviance_train.dropna())) + 1, np.array(deviance_train.dropna()['mean']),
         'b-',
         label='Training Set Deviance', )

plt.fill_between(np.arange(len(deviance_train.dropna())) + 1,
                 np.array(deviance_train.dropna()['mean']) + train_err,
                 np.array(deviance_train.dropna()['mean']) - train_err,
                 alpha=0.2)

plt.plot(np.arange(len(deviance_test.dropna())) + 1, np.array(deviance_test.dropna()['mean']),
         'r-',
         label='Test Set Deviance')

plt.fill_between(np.arange(len(deviance_test.dropna())) + 1,
                 np.array(deviance_test.dropna()['mean']) + test_err,
                 np.array(deviance_test.dropna()['mean']) - test_err,
                 alpha=0.2, color='red')

plt.legend(loc='upper right')
plt.xlabel('Boosting Iterations', fontsize=14)
plt.ylabel('Deviance', fontsize=14)
plt.savefig('./' + args.prefix + '.analysis/mean_deviance.pdf',
            dpi=300, bbox_inches='tight')
plt.close()

print("Plotting mean ROC curve and error region")
print()
cols = roc_fpr_all.columns
roc_fpr_all['mean'] = roc_fpr_all.apply(lambda x: np.mean(x[cols]), 1)
roc_fpr_all['std'] = roc_fpr_all.apply(lambda x: np.std(x[cols]), 1)

roc_fpr_all.to_csv('./' + args.prefix + '.analysis/classifier_metrics/roc_false_positive_rate.tsv',
                   sep='\t')

cols = roc_tpr_all.columns
roc_tpr_all['mean'] = roc_tpr_all.apply(lambda x: np.mean(x[cols]), 1)
roc_tpr_all['std'] = roc_tpr_all.apply(lambda x: np.std(x[cols]), 1)

roc_tpr_all.to_csv('./' + args.prefix + '.analysis/classifier_metrics/roc_true_positive_rate.tsv',
                   sep='\t')

cols = roc_auc_all.columns
roc_auc_all['mean'] = roc_auc_all.apply(lambda x: np.mean(x[cols]), 1)
roc_auc_all['std'] = roc_auc_all.apply(lambda x: np.std(x[cols]), 1)
roc_auc_all = roc_auc_all.reset_index()
roc_auc_all['index'] = 'ROC AUC'
roc_auc_all = roc_auc_all.set_index('index')

roc_auc_all.to_csv('./' + args.prefix + '.analysis/classifier_metrics/roc_area_under_curve.tsv',
                   sep='\t')

plt.figure(figsize=(6, 6))
plt.subplot(1, 1, 1)

fpr_err = np.array(roc_fpr_all['std'])
tpr_err = np.array(roc_tpr_all['std'])
m_auc = str(round(np.array(roc_auc_all['mean'])[0], 3))
m_auc_std = str(round(np.array(roc_auc_all['std'])[0], 3))

plt.plot(roc_fpr_all['mean'], roc_tpr_all['mean'], '-', color='black',
         label='Mean ROC curve\n(AUC=' + m_auc + ', STD=' + m_auc_std + ')')

# plt.fill_between(roc_fpr_all['mean'],
#                 roc_tpr_all['mean']+tpr_err, roc_tpr_all['mean']-tpr_err,
#                 color='blue', alpha=0.2)

# plt.fill_betweenx(roc_tpr_all['mean'],
#                  roc_fpr_all['mean']+fpr_err, roc_fpr_all['mean']-fpr_err,
#                  color='red', alpha=0.2)

plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate', fontsize=14)
plt.ylabel('True Positive Rate', fontsize=14)
plt.title('Mean receiver operating characteristic', fontsize=14)
plt.legend(loc="lower right")
plt.savefig('./' + args.prefix + '.analysis/classifier_metrics/mean_roc.pdf',
            dpi=300, bbox_inches='tight')
plt.close()

print("Plotting mean Precision-Recall curve and error region")
print()
cols = pre_all.columns
pre_all['mean'] = pre_all.apply(lambda x: np.mean(x[cols]), 1)
pre_all['std'] = pre_all.apply(lambda x: np.std(x[cols]), 1)

pre_all.to_csv('./' + args.prefix + '.analysis/classifier_metrics/precision.tsv', sep='\t')

cols = rec_all.columns
rec_all['mean'] = rec_all.apply(lambda x: np.mean(x[cols]), 1)
rec_all['std'] = rec_all.apply(lambda x: np.std(x[cols]), 1)

rec_all.to_csv('./' + args.prefix + '.analysis/classifier_metrics/recall.tsv', sep='\t')

cols = roc_auc_all.columns
av_pre_all['mean'] = av_pre_all.apply(lambda x: np.mean(x[cols]), 1)
av_pre_all['std'] = av_pre_all.apply(lambda x: np.std(x[cols]), 1)
av_pre_all = av_pre_all.reset_index()
av_pre_all['index'] = 'P-R AUC'
av_pre_all = av_pre_all.set_index('index')

av_pre_all.to_csv(
    './' + args.prefix + '.analysis/classifier_metrics/precision_recall_area_under_curve.tsv',
    sep='\t')

plt.figure(figsize=(6, 6))
plt.subplot(1, 1, 1)

pre_err = np.array(pre_all['std'])
rec_err = np.array(rec_all['std'])
m_auc = str(round(np.array(av_pre_all['mean'])[0], 3))
m_auc_std = str(round(np.array(av_pre_all['std'])[0], 3))

plt.plot(rec_all['mean'], pre_all['mean'], '-', color='black',
         label='Mean Precision-Recall curve\n(AUC=' + m_auc + ', STD=' + m_auc_std + ')')

# plt.fill_between(rec_all['mean'],
#                 pre_all['mean']+pre_err, pre_all['mean']-pre_err,
#                 color='blue', alpha=0.2)

# plt.fill_betweenx(pre_all['mean'],
#                  rec_all['mean']+rec_err, rec_all['mean']-rec_err,
#                  color='red', alpha=0.2)

plt.xlim([-0.05, 1.05])
plt.ylim(ymax=1.05)
plt.xlabel('Recall', fontsize=14)
plt.ylabel('Precision', fontsize=14)
plt.title('Precision-Recall', fontsize=14)
plt.legend(loc="lower left")
plt.savefig(
    './' + args.prefix + '.analysis/classifier_metrics/mean_precision_recall.pdf',
    dpi=300, bbox_inches='tight')
plt.close()

print(
    "Plotting barcharts for mean classifier scores, confusion matrix values and associated metrics")
print()

df_list = [conf_df, msr, mtc]
title_list = ['Identified classes from confusion matrix',
              'MSE and aMI scores from classifier',
              'Metrics estimated from confusion matrix']
save_list = ['mean_confusion_matrix_classes',
             'mse_ami_scores_from_classifier',
             'classifier_metrics_from_confusion_matrix']

for run_i in range(len(df_list)):
    plot_barcharts(df_list[run_i], title_list[run_i], save_list[run_i])
    plt.savefig('./' + args.prefix + '.analysis/classifier_metrics/' + save_list[
        run_i] + '.pdf', dpi=300, bbox_inches='tight')
    plt.close()

mtc_ms = get_mean_and_std(mtc)
msr_ms = get_mean_and_std(msr)
    
df_cat = pd.concat([mtc_ms[['mean','std']],
                    msr_ms[['mean','std']],
                    roc_auc_all[['mean','std']],
                    av_pre_all[['mean','std']]])

df_cat.to_csv('./' + args.prefix + '.analysis/overall_classifier_metrics.tsv', sep='\t')
df_cat.to_excel('./' + args.prefix + '.analysis/overall_classifier_metrics.xlsx')

N = df_cat.shape[0]
means = df_cat['mean'].round(2)
std = df_cat['std']

ind = np.arange(N)  # the x locations for the groups
width = 0.75  # the width of the bars

sns.set_style('ticks')

fig, ax = plt.subplots()
rects1 = ax.bar(ind, means, width, 
                edgecolor='black', linewidth=0.5, color='white', 
                yerr=std, ecolor='black')

# add some text for labels, title and axes ticks
ax.set_ylabel('Value', fontsize=14)
ax.set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0], fontsize=12)
ax.set_xticks(ind)
ax.set_xticklabels(df_cat.index, rotation=45, fontsize=12)
ax.grid(False)

rects = ax.patches

# Now make some labels
labels = ["%.2f" % i for i in means]

for rect, label in zip(rects, labels):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2, height+0.01, label, ha='center', va='bottom')
    

sns.despine(top=True, right=True, left=False, bottom=False, offset=None, trim=False)
    
plt.savefig(
    './' + args.prefix + '.analysis/overall_classifier_metrics.pdf',
    dpi=300, bbox_inches='tight')

print ("Plotting overall feature importance and KS-statistics data")
print()

sns.set_style('ticks')

if importance.shape[0] >= 10:
    N = 10
else:
    N = importance.shape[0]

importance_ind = importance.sort_values(by=['mean_rel_importance'], ascending=False).set_index("Feature").head(N)

ind = np.arange(N)  # the x locations for the groups
width = 0.25  

rel = importance_ind['mean_rel_importance'].head(N)
rel_err = importance_ind['std_rel_importance'].head(N)
names = (str(rel.index)).split("/")[-1]

#    importance_ind['-log10(q-value)'] = importance_ind['adj_pval_0_vs_1'].apply(lambda x: (np.log10(x))*-1, 1)
ks = importance_ind['ks'].head(N)

fig = plt.figure(figsize=(6,6)) # Create matplotlib figure
ax = fig.add_subplot(111) # Create matplotlib axes
ax2 = ax.twiny() # Create another axes that shares the same x-axis as ax.

rel.plot(kind='barh', edgecolor='black', linewidth=0.5, color='white', 
         label='Relative importance', ax=ax2, width=width, position=1,
         xerr=rel_err, ecolor='black')

ks.plot(kind='barh', edgecolor='black', linewidth=0.5, color='lightgrey', 
        label='KS test', ax=ax, width=width, position=0)

ax2.legend(loc='center left', bbox_to_anchor=(0.05, -0.175), fontsize=12)
ax2.tick_params(axis='both', which='major', labelsize=14)
ax2.grid(b=False)
ax2.patch.set_visible(False)
ax2.set_xlabel('Relative importance', size=14)
ax2.tick_params(which='major', axis='x', pad=5)
ax2.xaxis.labelpad = 10

ax.legend(loc='center left', bbox_to_anchor=(0.55, -0.175), fontsize=12)
ax.set_xlabel('KS test', size=14)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.grid(b=False)
ax.patch.set_visible(False)

ax2.tick_params(top=False, bottom=False, left=False, right=False, labelleft=True, labeltop=True)
ax.tick_params(top=False, bottom=False, left=False, right=False, labelleft=True, labelbottom=True)

plt.gca().invert_yaxis()
plt.tick_params(axis='both', which='major', labelsize=14)

sns.despine(offset=5, trim=True, ax=ax)
sns.despine(offset=5, trim=True, ax=ax2)

plt.savefig(
    './' + args.prefix + '.analysis/classifier_results.pdf',
    dpi=300, bbox_inches='tight')
    
print("Zipping folders from every run in a single file")
print()

Popen(
    'tar -zcf ./' + args.prefix + '.analysis/classifier_metrics/run_files.tar.gz ./' + args.prefix + '.analysis/classifier_metrics/run_*/',
    shell=True)

print("Analysis complete. Thanks for using biofeatures.")
print()
sys.exit()
