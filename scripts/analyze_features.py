#!/usr/python

#### Load the required packages

import sys
import argparse
from argparse import ArgumentParser
import os
import time
import subprocess
from subprocess import PIPE,call,Popen
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
matplotlib.use('Agg')
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
import multiprocessing as mp
pd.options.mode.chained_assignment = None  # default='warn'


plt.ioff()

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

parser.add_argument('-ts',"--twoSamp",dest="twoSamp",
                  action="store_true", default=False,
                  help="Use this flag if you want to compare two sets of samples against each other, requires the usage of -b2. In this case, the second bed file will be considered as 'negative' (or background) for the classification step. Default: False")

parser.add_argument('-b', '--bed', dest="bed_file", 
                  help="BED file with exons/regions of interest.", required=True)

parser.add_argument('-b2', '--bed2', dest="bed_file_2", 
                  help="Second BED file with exons/regions of interest.", required=False)

parser.add_argument('-ref', '--refference_bed', dest="refference_bed", 
                  help="BED file with annotated exons created by 'build_datamatrix.py'", required=True)

parser.add_argument('-m', '--matrix', dest="matrix", 
                  help="Data matrix with biological features created by 'build_datamatrix.py'", required=True)

parser.add_argument('-o', '--outfile', dest="prefix", #type=str,
                  help="prefix for use on the output files", metavar="prefix", required=True)

parser.add_argument('-filter', '--filter_columns', dest="filter_out", default=False,
                  help="Text file containing a comma-separated list with names of the columns to be removed from the dataframe in the analysis. Default: False", metavar="filter_out.txt", required=False)

parser.add_argument('-select', '--select_columns', dest="filter_in", default=False,
                  help="Text file containing a comma-separated list with names of the columns in dataframe to be used in the analysis. Default: False", metavar="filter_in.txt", required=False)

parser.add_argument("-padj", '--p_adjust', dest="padj", default='bonferroni', help="Type of p-value correction used after Kolmogorov-Smirnov test, available options are: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' and 'none'. Default:'bonferroni'", type=str, metavar='padj', required=False)

parser.add_argument('-F', '--F_bedtools', dest="F_bedtools", default=1, help="-F options passed to BedTools for intersection of reference BED (-ref, A) with input BED (-b, B). Default: 1", metavar='F', required=False)

parser.add_argument('-f','--f_bedtools', dest="f_bedtools", default=1, help="-f options passed to BedTools for intersection of reference BED (-ref, A) with input BED (-b, B). Default: 1",  metavar='f',required=False)

parser.add_argument("--skipStat",dest="skip_stat",
                  action="store_true", default=False,
                  help="Use this flag if you want to skip the Kolmogorov-Smirnov statistical analysis and cumulative distribution plotting step. Default: False")

parser.add_argument("--no-plotCDF",dest="dont_plot_cdf",
                  action="store_true", default=False,
                  help="Use this flag if you want to skip plotting CDF graphs for each feature in the matrix. Default: False")

parser.add_argument('--sig-only-CLF', dest="ks_filter", action="store_true", default=False, help="Use only the statistically significant features (found by KS test) in the classification step. Useful for filtering large data matrices to reduce noise error and increase accuracy, however it may require some manual parameter tunning (with -params 'file' option). Can use the '-pth' option to select the threshold of significante for feature selection. Default: False", required=False)

parser.add_argument("-pth", '--p_adjust_threshold', dest="p_th", default=0.05, help="Threshold of adjusted p-value for significance. Only significantly different features are passed down to the classifier for group separation. Default: 0.05", type=float, required=False)

parser.add_argument("-runs",'--number_of_runs', dest="runs", default=(500), help="Number of times (repetitions) to run the classification step. Default:500", type=int, metavar='INT')

parser.add_argument("-nsample",'--random_sample_size', dest="nsample", type=int, metavar='INT', default=(1), help="Relative size of randomly sampled exons in comparisson to input exons. Default:1 (i.e. 1x the amount of input exons)")

parser.add_argument("-tsize",'--train_size', dest="train_size", default=(0.80), help="Fraction of sample used for training the classifier model. The remaining sample pool will be used for testing the classifier. Default: 0.80", type=float)

parser.add_argument("-params", '--gbcl_parameters', dest="clf_params", default='preset', help="Type of parameter selection to be used by the classifier. Available options are: 'optimize' (runs an optimization step with GridSearchCV before every run), 'default' (uses default parameters from GradientBoostClassifier), 'preset' (uses the preset parameters which were used for the analysis of the human dataset shown in the article) and 'file' (take in input txt file with a dictionary-like struture with classifier parameters, requires the use of -pf option). Default:'preset'", type=str, metavar='params', required=False)

parser.add_argument("-pf", '--param_file', dest="param_file", help="Input text with with dictionary-like structure with parameter options for GradientBoostClassifier. Ex. {'n_estimators':300,'loss':'deviance',...}", metavar='file', required=False)

parser.add_argument("--noCLF",dest="dont_run_clf", 
                  action="store_true", default=False,
                  help="Use this flag if you want to skip the classifying with GradientBoost. Default: False")

parser.add_argument("--ncores", dest="ncores", default=(mp.cpu_count()-1), help="Number of CPU cores used to multiple jobs on the classifier. Default:(ALL_CORES)-1", type=int, metavar='INT')

args = parser.parse_args()

##Define the functions which will be used during the analysis

def group_matrices_one_sample(bt, bt_a, matrix, F=args.F_bedtools, f=args.f_bedtools, s=True):
    
    int_a = bt.intersect(bt_a, 
                         s=s, 
                         F=F, 
                         f=f, 
                         sorted=True).to_dataframe().drop_duplicates()
    
    matrix['group'] = 0
    matrix['group'][matrix['name'].isin(int_a['name'])] = 1
    
    return matrix

def group_matrices_two_samples(bt, bt_a, bt_b, matrix, F=args.F_bedtools, f=args.f_bedtools, s=True):
    
    int_a = bt.intersect(bt_a, 
                         s=s, 
                         F=F, 
                         f=f, 
                         sorted=True).to_dataframe().drop_duplicates()
    
    int_b = bt.intersect(bt_b, 
                         s=s, 
                         F=F, 
                         f=f, 
                         sorted=True).to_dataframe().drop_duplicates()
    
    matrix['group'] = 0
    matrix['group'][(matrix['name'].isin(int_a['name'])) & 
                    ~(matrix['name'].isin(int_b['name']))] = 1
    
    matrix['group'][(matrix['name'].isin(int_b['name'])) & 
                    ~(matrix['name'].isin(int_a['name']))] = 2
    
    matrix['group'][(matrix['name'].isin(int_b['name'])) & 
                    (matrix['name'].isin(int_a['name']))] = 3
    
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

def plot_barcharts(df, title, save):
    df_t = df.T
    cols = df_t.columns

    df_t['mean'] = df_t.apply(lambda x: np.mean(x[cols]),1)
    df_t['std'] = df_t.apply(lambda x: np.std(x[cols]),1)

    df_t.to_csv('./'+args.prefix+'_results/'+save+'.tsv', sep='\t')
    
    N = df_t.shape[0]
    means = df_t['mean']
    std = df_t['std']

    ind = np.arange(N)  # the x locations for the groups
    width = 0.5       # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, means, width, color='lightgrey', yerr=std, ecolor='black')

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Value', fontsize=14)
    ax.set_title(title, fontsize=14)
    ax.set_xticks(ind+width/2)
    ax.set_xticklabels(df_t.index, rotation='vertical')
    
def plot_barchart_importance(df):
    importance_ind = df.set_index("Feature")

    N = 10
    ind = np.arange(N)  # the x locations for the groups
    width = 0.5    

    rel = importance_ind.sort_values(by=['mean_rel_importance'], ascending=False)['mean_rel_importance'].head(N)
    rel_err = importance_ind.sort_values(by=['mean_rel_importance'], ascending=False)['std_rel_importance'].head(N)

    fig, ax = plt.subplots()
    ax.bar(ind, rel, width, color='lightgrey', yerr=rel_err, ecolor='black')
    ax.set_ylabel('Relative importance', fontsize=14)
    ax.set_title('Mean relative importance values', fontsize=14)
    ax.set_xticks(ind+width/2)
    ax.set_xticklabels(rel.index, rotation='vertical', fontsize=14)
    plt.savefig('./'+args.prefix+'_results/classifier_plots/mean_relative_importance.pdf',dpi=300, bbox_inches='tight')
    plt.close()

    raw = importance_ind.sort_values(by=['mean_raw_importance'], ascending=False)['mean_raw_importance'].head(N)
    raw_err = importance_ind.sort_values(by=['mean_raw_importance'], ascending=False)['std_raw_importance'].head(N)

    fig, ax = plt.subplots()
    ax.bar(ind, raw, width, color='lightgrey', yerr=raw_err, ecolor='black')
    ax.set_ylabel('Importance', fontsize=14)
    ax.set_title('Mean importance values', fontsize=14)
    ax.set_xticks(ind+width/2)
    ax.set_xticklabels(raw.index, rotation='vertical', fontsize=14)
    plt.savefig('./'+args.prefix+'_results/classifier_plots/mean_importance.pdf',dpi=300, bbox_inches='tight')
    plt.close()
    
##Load the datamatrix generatade by "buildadatamatrix.py"

print
print "Loading datamatrix"
print

matrix = pd.concat(pd.read_table(args.matrix,iterator=True, chunksize=10000), 
                          ignore_index=True).set_index('name').drop_duplicates()#.reset_index()

#Filter in/out columns in the dataframe

if args.filter_out == False:
    pass
else:
    print "Filtering out columns"
    print
    out_cols = open(str(args.filter_out)).read().split(',')
    out_cols = [w.replace('\n', '') for w in out_cols]
    matrix = matrix.drop(out_cols,1)
    
if args.filter_in == False:
    pass
else:
    print "Selecting columns"
    print
    in_cols = open(str(args.filter_in)).read().split(',')
    in_cols = [w.replace('\n', '') for w in in_cols]
    matrix = matrix[in_cols]

##Load the bed file created with all exons

matrix = matrix.reset_index()

print "Loading annotated exons in refference file"
print 

all_exons = BedTool(args.refference_bed).sort()

##Load the dataset of interest to be analyzed (needs to be a SORTED BED)

print "Load selected exons for biological feature analysis"
print

if args.twoSamp == False:
    pass
elif args.twoSamp == True:
    print "Running in two sample mode"
    print
    print "Loading bed file for sample 1"
    print

sel_exons = BedTool(args.bed_file).sort()

if args.twoSamp == True:
    print "Loading bed file for sample 2"
    print
    sel_exons_2 = BedTool(args.bed_file_2).sort()
else:
    pass

##Intersect the exons found in the analysis to get groups 1 (positive) and 0 (negative) in the matrix

print "Finding input exons in the matrix and selecting groups"
print

if args.twoSamp == False:
    matrix = group_matrices_one_sample(all_exons, sel_exons, matrix).set_index('name')
elif args.twoSamp == True:
    matrix = group_matrices_two_samples(all_exons, sel_exons, sel_exons_2, matrix).set_index('name')
    

##Run statistical analysis on the dataset with the "get_statistical_data_for_features" function:

Popen('mkdir -p ./'+args.prefix+'_results', shell=True)

if args.skip_stat == False:
    print "Starting statistical analysis"
    print

    st = get_statistical_data_for_features(matrix, correction_type=args.padj)
    st.to_csv('./'+args.prefix+'_results/statistical_analysis_output.tsv', sep='\t', index=False)
    
    print "Finished statistical analysis"
    print
    ##########################################
    ## Plot CDF for the features
    if args.dont_plot_cdf == False:
        print "Output CDF plots for each features in matrix"
        print 
        Popen('mkdir -p ./'+args.prefix+'_results/CDF_plots', shell=True)

        for i in range(len(st)):
            feature = st['Feature'][i]
            name = (st['Feature'][i]).split("/")[-1]
            plt.figure(figsize=(8,8))              
            plot_cdf(matrix[matrix['group'] == 0][str(feature)].values, bins=100, 
                     label='Exons in 0', c='black', linewidth=1.5, linestyle='solid');
            plot_cdf(matrix[matrix['group'] == 1][str(feature)].values, bins=100, 
                     label='Exons in 1', c='black', linewidth=1.5, linestyle='dashed');

            if args.twoSamp == False:
                pass
            elif args.twoSamp == True:
                plot_cdf(matrix[matrix['group'] == 2][str(feature)].values, bins=100, 
                         label='Exons in 2', c='black', linewidth=1.5, linestyle='dotted');

                plot_cdf(matrix[matrix['group'] == 3][str(feature)].values, bins=100, 
                         label='Exons in both 1 and 2', c='black', linewidth=1.5, linestyle='-.');

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

            plt.savefig('./'+args.prefix+'_results/CDF_plots/'+name+'.pdf', dpi=300, bbox_inches='tight')
            plt.close()

        print "Finished CDF plots for features in matrix"
        print
    elif args.dont_plot_cdf == True:
        pass
    
    if args.ks_filter == False:
        pass
    elif args.ks_filter == True:
        if args.twoSamp == False:
            print "Filtering statistically significantly features for classification step"
            print
            sig_only = st[st['adj_pval_0_vs_1'] <= float(args.p_th)]['Feature'].tolist()
            sig_only.append('group')
            matrix = matrix[sig_only]
        elif args.twoSamp == True:
            print "Filtering statistically significantly features for classification step"
            print
            sig_only = st[st['adj_pval_1_vs_2'] <= float(args.p_th)]['Feature'].tolist()
            sig_only.append('group')
            matrix = matrix[sig_only]
    
elif args.skip_stat == True:
    pass
    
if args.dont_run_clf == True:
    print "Analysis complete. Thanks for using BioFeatureFinder."
    print
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

Popen('mkdir -p ./'+args.prefix+'_results/classifier_plots', shell=True)

if args.skip_stat == False:
    importance = st.copy()
elif args.skip_stat == True:
    importance = pd.DataFrame()
    importance['Feature'] = matrix.drop('group',1).columns
    
deviance_train = pd.DataFrame()
deviance_test = pd.DataFrame()

msr = pd.DataFrame(columns=['aMI','MSE'])

mtc = pd.DataFrame(columns=['Accuracy','P.P.V.','N.P.V.','Sensitivity','Specificity'])

conf_df = pd.DataFrame(columns=['T.N.','T.P.','F.N.','F.P.'])

roc_tpr_all = pd.DataFrame()
roc_fpr_all = pd.DataFrame()
roc_auc_all = pd.DataFrame()

pre_all = pd.DataFrame()
rec_all = pd.DataFrame()
av_pre_all = pd.DataFrame()

if args.twoSamp == True:
    matrix = matrix[(matrix['group'] == 2) | (matrix['group'] == 1)]
    matrix['group'] = matrix['group'].apply(lambda x: int(str(x).replace('2','0')))
else:
    pass

for i in range(len(runs)):
    run_id = str(runs[i]+1)
    if args.twoSamp == True:
        df_cl = matrix[matrix['group'] != 0]
        df_zero = matrix[matrix['group']==0]#.sample(df_cl.shape[0]*4)
        Z = pd.concat([df_cl,df_zero]).drop_duplicates()
    else:
        df_cl = matrix[matrix['group'] != 0]
        df_zero = matrix[matrix['group']==0].sample(df_cl.shape[0]*sample_size)
        Z = pd.concat([df_cl,df_zero]).drop_duplicates()
    
    Popen('mkdir -p ./'+args.prefix+'_results/classifier_plots/run_'+run_id, shell=True)
    
    print "Run ID: "+str(run_id)+', Input size: '+str(df_cl.shape[0])+', Background size: '+str(df_zero.shape[0])
    print 
    print "Starting classification analysis with GradientBoost. This may take a long time depending on the size of the dataset"
    print

    ##The first step is separating the groups from the data in different arrays

    print "Loading data in arrays"
    print

    y = np.array(Z.group).astype('int')
    X = np.array(Z.drop('group', 1))

    ##Then, we split the total data in 2 groups: training group (which we will train our model and optimize the parameters) 
    ## and the test group, which we will use for evaluating the performance of the classifier. The default test size is 0.25.

    print "Splitting train and test datasets"
    print

    X_train, X_test, y_train, y_test = train_test_split(X,
                                                        y,
                                                        train_size=train_size)

    ##For parameter selection to be passed to the classifier, there are 4 options: 
    ## 1. "optimize" - Run the optimization step which performs a GridSearchCV on several combinations of parameters 
    ## to find the best performing one. 
    ## 2. "default" - Run the default classfier parameters 
    ## 3. "preset" - Use the preset parameters that was used in the human dataset testing shown in the article. 
    ## 4. Input a text file containing a dictionary-like list of parameters of your own choice.
    if gbclf_params == 'optimize':

        ##Since we are using the GradientBooost Classfier, the first step is to find how many features are 
        ##optimal for group classification. To do that, we use recursive feature elimination with cross-validation.
        ##For this step, we will use the StratifiedKFold approach with 10-fold cross validation. The "accuracy" 
        ##scoring is proportional to the number of correct classifications

        clf = GradientBoostingClassifier(warm_start=True)

        print "Optimizing number of features with recursive elimination using stratified 10-fold cross-validation loop"
        print

        rfecv = RFECV(estimator=clf, step=1, cv=StratifiedKFold(10),
                      scoring='accuracy', n_jobs=args.ncores)

        rfecv.fit(X, y)

        print("Optimal number of features : %d" % rfecv.n_features_)
        print

        best_n_features = rfecv.n_features_

        # Plot number of features VS. cross-validation scores ('1st metric')
        print "Plotting N.features x Cross-validation scores curve"
        print

        plt.figure()
        plt.xlabel("Number of features selected")
        plt.ylabel("Cross validation score")
        plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_, lw=2)
        plt.savefig('./'+args.prefix+'_results/classifier_plots/run_'+run_id+'/run_'+run_id+'_feature_recursive_elimination.pdf', dpi=300, bbox_inches='tight')
        plt.close()

        print "Running GridSearchCV to optimize remaining parmeters"
        print 

        grid = {'n_estimators':[500,1000,2000],
                'max_depth':[6,8,10],
                'min_samples_split':[0.01],
                'min_samples_leaf':[0.001],
                'max_features':[best_n_features],
                'learning_rate':[0.01,0.05], 
                'loss':['deviance'],
                'subsample':[0.8,0.6],
                'random_state': [1]}

        gclf = GridSearchCV(estimator=clf, param_grid=grid, n_jobs=args.ncores, cv=10)

        gclf.fit(X_train,y_train)
        bp = gclf.best_params_
        
        with open ('./'+args.prefix+'_results/classifier_plots/run_'+run_id+'/run_'+run_id+'_best_params.txt', 'w') as fp:
            fp.write('{')
            for p in bp.items():
                fp.write("'%s':%s,\n" % p)
            fp.write('}')

        print "Confusion matrix: "
        print confusion_matrix(y_test, gclf.best_estimator_.predict(X_test))
        print
        print "GridSeachCV best score:"
        print gclf.best_score_
        print
        print "GridSearchCV best params"
        print gclf.best_params_
        print
    if gbclf_params == 'default':
        print "Using default parameters."
        print
    if gbclf_params == 'preset':
        print "Using preset parameters."
        print
        bp = {'learning_rate': 0.01,
              'loss': 'deviance',
              'max_depth': 8,
              'max_features': 'sqrt',
              'min_samples_leaf': 0.001,
              'min_samples_split': 0.01,
              'n_estimators': 1000,
              'random_state': 1,
              'subsample': 0.8}
    if gbclf_params == 'file':
        print "Using input file as param dictionary"
        print
        bp = eval(open(str(args.param_file)).read())

    print "Fitting the GradientBoost model with the training set"
    print
    
    names = df_cl.drop('group',1).columns

    if gbclf_params == 'default':
        clf = GradientBoostingClassifier(warm_start=True)
    else:
        clf = GradientBoostingClassifier(warm_start=True,**bp)

    clf.fit(X_train, y_train)
    
    print "Extracting feature importance for run "+run_id+" and merging with statistical data"
    print
    
    a = pd.DataFrame(zip(clf.feature_importances_, 
        np.argsort(np.argsort(clf.feature_importances_)), names)).rename(columns={0:'run_'+run_id+'_raw_importance',
                                                                                  1:'Index',
                                                                                  2:'Feature'}
                                                                        ).sort_values('Index', ascending=False)
    
    max_importance = a['run_'+run_id+'_raw_importance'].max()
    a['run_'+run_id+'_rel_importance'] = 100.0 * (a['run_'+run_id+'_raw_importance'] / max_importance)
    importance = importance.merge(a.drop('Index',1), on='Feature')
    
    print "Plotting partial dependance and deviance scores from classifier"
    print
       
    #Select top 5 most important features and plot partial dependance and interaction plots
    
    do_these = a.sort_values(by='Index', ascending=False).index[:4].values

    features = list(do_these)
    features.append((do_these[0], do_these[1]))
    fig, axs = plot_partial_dependence(clf, X_train, features, feature_names=names,
                                       n_jobs=1, grid_resolution=15, figsize=(14,8),
                                       label=0)

    plt.subplots_adjust(top=0.9)  # tight_layout causes overlap with suptitle

    plt.tight_layout()
    plt.savefig('./'+args.prefix+'_results/classifier_plots/run_'+run_id+'/run_'+run_id+'_partial_dependance.pdf', dpi=600,
               bbox_inches='tight')
    plt.close()
    
    ##Plot deviance x number os iteration
    ###############################################################################
    params = clf.get_params()
    test_score = np.zeros((params['n_estimators'],), dtype=np.float64)

    for i, y_pred in enumerate(clf.staged_decision_function(X_test)):
        test_score[i] = clf.loss_(y_test, y_pred)

    plt.figure(figsize=(6, 6))
    plt.subplot(1, 1, 1)
    plt.title('Deviance')
    plt.plot(np.arange(params['n_estimators']) + 1, clf.train_score_, 'b-',
             label='Training Set Deviance')
    plt.plot(np.arange(params['n_estimators']) + 1, test_score, 'r-',
             label='Test Set Deviance')
    plt.legend(loc='upper right')
    plt.xlabel('Boosting Iterations')
    plt.ylabel('Deviance')
    plt.savefig(''+args.prefix+'_results/classifier_plots/run_'+run_id+'/run_'+run_id+'_deviance.pdf', dpi=600,
                bbox_inches='tight')
    plt.close()
    
    deviance_train['run_'+run_id] = clf.train_score_
    deviance_test['run_'+run_id] = test_score
    
    print "Run trained classifier in test subset from data"
    print
    predicted_values = clf.predict(X_test)
    mse = mean_squared_error(y_test, predicted_values)
    ami = adjusted_mutual_info_score(y_test, predicted_values)
    conf = confusion_matrix(y_test, predicted_values)
    spe = round((float(conf[0,0])/(conf[1,0]+conf[0,0]))*100,2)
    sen = round((float(conf[1,1])/(conf[0,1]+conf[1,1]))*100,2)
    npv = round((float(conf[0,0])/(conf[0,0]+conf[0,1]))*100,2)
    ppv = round((float(conf[1,1])/(conf[1,0]+conf[1,1]))*100,2)
    acc = round((float(conf[0,0]+conf[1,1])/(conf[0,0]+conf[1,0]+conf[0,1]+conf[1,1]))*100,2)
    
    msr_r = pd.DataFrame(data=[[ami,mse]], columns=['aMI','MSE'])    
    
    msr = msr.append(msr_r).reset_index().drop('index',1)
    
    mtc_r = pd.DataFrame(data=[[acc,ppv,npv,sen,spe]], 
                         columns=['Accuracy','P.P.V.','N.P.V.','Sensitivity','Specificity'])
    
    mtc = mtc.append(mtc_r).reset_index().drop('index',1)
    
    conf_df_r = pd.DataFrame(data=[[conf[0,0],conf[1,1],conf[1,0],conf[0,1]]],
                             columns=['T.N.','T.P.','F.N.','F.P.'])
    
    conf_df = conf_df.append(conf_df_r).astype(int).reset_index().drop('index',1)

    print "Classifier scores"
    print("Mean Squared Error (MSE): %.4f" % mse)
    print("adjusted Mutual Information (aMI): %.4f" % ami)
    print 
    print "Confusion matrix"
    print(conf)
    print
    print 'Accuracy: '+str(acc)+'%'
    print 'Positive predictive value (Precision): '+str(ppv)+'%'
    print 'Negative predictive value: '+str(npv)+'%'
    print 'Sensitivity (Recall): '+str(sen)+'%'
    print 'Specificity: '+str(spe)+'%'
    print
    
    ##For creating ROC and precision curves, we need to binarize the output and get score functions
    print "Binarizing output"
    print
    y_class = label_binarize(y_test, classes=[0,1])
    n_classes = y_class.shape[1]
    y_score = clf.fit(X_train, y_train).decision_function(X_test)
    
    ###Now, let's calculate the ROC score for the classification
    
    print "Computing ROC score for classifier"
    print
    
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    # Compute micro-average ROC curve and ROC area
    fpr["micro_run_"+run_id], tpr["micro_run_"+run_id], _ = roc_curve(y_test.ravel(), y_score.ravel())
    roc_auc["micro_run_"+run_id] = auc(fpr["micro_run_"+run_id], tpr["micro_run_"+run_id])
      
    fpr_r = pd.DataFrame(data=fpr["micro_run_"+run_id], columns=['F.P.R.'])
    fpr_r.to_csv('./'+args.prefix+'_results/classifier_plots/run_'+run_id+'/run_'+run_id+'_roc_fpr.tsv', sep='\t')
    
    tpr_r = pd.DataFrame(data=tpr["micro_run_"+run_id], columns=['T.P.R.'])
    tpr_r.to_csv('./'+args.prefix+'_results/classifier_plots/run_'+run_id+'/run_'+run_id+'_roc_tpr.tsv', sep='\t')
        
    roc_fpr_all = pd.concat([roc_fpr_all,fpr_r], ignore_index=True, axis=1)#.fillna(1)
    roc_tpr_all = pd.concat([roc_tpr_all,tpr_r], ignore_index=True, axis=1)#.fillna(1)
    roc_auc_all['run_'+run_id] = np.array([roc_auc['micro_run_'+run_id]])

    # Plot ROC curve

    plt.plot(fpr["micro_run_"+run_id], tpr["micro_run_"+run_id],
            label='ROC curve (AUC = {0:0.2f})'
                  ''.format(roc_auc["micro_run_"+run_id]), linewidth=2)

    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate', size=20)
    plt.ylabel('True Positive Rate', size=20)
    plt.title('Receiver operating characteristic', size=20)
    plt.legend(loc="lower right", fontsize=18)
    plt.savefig('./'+args.prefix+'_results/classifier_plots/run_'+run_id+'/run_'+run_id+'_roc_curve.pdf', dpi=300,
                bbox_inches='tight')
    plt.close()
    
    print "Plotting Precision-Recall curve"
    print
    # Compute Precision-Recall and plot curve
    precision = dict()
    recall = dict()
    average_precision = dict()
    
    precision["micro_run_"+run_id], recall["micro_run_"+run_id], _ = precision_recall_curve(y_test,y_score)
    average_precision["micro_run_"+run_id] = average_precision_score(y_test, y_score)
    
    pre_r = pd.DataFrame(data=precision["micro_run_"+run_id], columns=['Precision'])
    pre_r.to_csv('./'+args.prefix+'_results/classifier_plots/run_'+run_id+'/run_'+run_id+'_precision.tsv', sep='\t')
    
    rec_r = pd.DataFrame(data=recall["micro_run_"+run_id], columns=['Recall'])
    rec_r.to_csv('./'+args.prefix+'_results/classifier_plots/run_'+run_id+'/run_'+run_id+'_recall.tsv', sep='\t')
        
    pre_all = pd.concat([pre_all,pre_r], ignore_index=True, axis=1)#.fillna(0)
    rec_all = pd.concat([rec_all,rec_r], ignore_index=True, axis=1)#.fillna(0)
    av_pre_all['run_'+run_id] = np.array([average_precision['micro_run_'+run_id]])
    
    # Plot Precision-Recall curve

    plt.figure(figsize=(6,6))

    plt.clf()
    plt.plot(recall["micro_run_"+run_id], precision["micro_run_"+run_id], lw=2, color='navy',
             label='Precision-Recall curve (AUC={0:0.2f})'.format(average_precision["micro_run_"+run_id]))
    plt.title('Precision-Recall', size=20)
    plt.xlabel('Recall', size=20)
    plt.ylabel('Precision', size=20)
    plt.xlim([-0.05, 1.05])
    #plt.ylim([-0.05, 1.05])
    plt.ylim(ymax=1.05)
    plt.legend(loc="lower left", fontsize=18)
    plt.savefig('./'+args.prefix+'_results/classifier_plots/run_'+run_id+'/run_'+run_id+'_precision_recall_curve.pdf', dpi=300,
                bbox_inches='tight')
    plt.close()
    
    print '======================='
    print

print "Done with classification steps"
print
print "Calculating median importance values and standard deviation"
print 

cols = importance.filter(like='_raw_importance',axis=1).columns
importance['mean_raw_importance'] = importance.apply(lambda x: np.mean(x[cols]),1)
importance['std_raw_importance'] = importance.apply(lambda x: np.std(x[cols]),1)

cols = importance.filter(like='_rel_importance',axis=1).columns
importance['mean_rel_importance'] = importance.apply(lambda x: np.mean(x[cols]),1)
importance['std_rel_importance'] = importance.apply(lambda x: np.std(x[cols]),1)

importance.to_csv('./'+args.prefix+'_results/classifier_importance.tsv', sep='\t')

print "Plotting raw importance and relative importance barcharts"
print

plot_barchart_importance(importance)

print "Plotting mean deviance curve and error region"
print
cols = deviance_test.columns
deviance_test['mean'] = deviance_test.apply(lambda x: np.mean(x[cols]),1)
deviance_test['std'] = deviance_test.apply(lambda x: np.std(x[cols]),1)

deviance_test.to_csv('./'+args.prefix+'_results/classifier_test_deviance.tsv', sep='\t')

cols = deviance_train.columns
deviance_train['mean'] = deviance_train.apply(lambda x: np.mean(x[cols]),1)
deviance_train['std'] = deviance_train.apply(lambda x: np.std(x[cols]),1)

deviance_train.to_csv('./'+args.prefix+'_results/classifier_train_deviance.tsv', sep='\t')

plt.figure(figsize=(6, 6))
plt.subplot(1, 1, 1)
plt.title('Deviance')

train_err = np.array(deviance_train['std'])
test_err = np.array(deviance_test['std'])

plt.plot(np.arange(len(deviance_train)) + 1, np.array(deviance_train['mean']), 'b-',
         label='Training Set Deviance', )

plt.fill_between(np.arange(len(deviance_train)) + 1,
                 np.array(deviance_train['mean'])+train_err, 
                 np.array(deviance_train['mean'])-train_err,
                alpha=0.2)

plt.plot(np.arange(len(deviance_test)) + 1, np.array(deviance_test['mean']), 'r-',
         label='Test Set Deviance')

plt.fill_between(np.arange(len(deviance_test)) + 1,
                 np.array(deviance_test['mean'])+test_err, 
                 np.array(deviance_test['mean'])-test_err,
                alpha=0.2, color='red')

plt.legend(loc='upper right')
plt.xlabel('Boosting Iterations', fontsize=14)
plt.ylabel('Deviance', fontsize=14)
plt.savefig('./'+args.prefix+'_results/classifier_plots/mean_deviance.pdf',dpi=300, bbox_inches = 'tight')
plt.close()

print "Plotting barcharts for mean classifier scores, confusion matrix values and associated metrics"
print

df_list = [conf_df,msr,mtc]
title_list = ['Identified classes from confusion matrix', 
              'MSE and aMI scores from classifier',
              'Metrics estimated from confusion matrix']
save_list = ['mean_confusion_matri_classes',
             'mse_ami_scores_from_classifier',
             'classifier_metrics_from_confusion_matrix']

for i in range(len(df_list)):
    plot_barcharts(df_list[i], title_list[i], save_list[i])
    plt.savefig('./'+args.prefix+'_results/classifier_plots/'+save_list[i]+'.pdf',dpi=300, bbox_inches='tight')
    plt.close()

    
print "Plotting mean ROC curve and error region"
print
cols = roc_fpr_all.columns
roc_fpr_all['mean'] = roc_fpr_all.apply(lambda x: np.mean(x[cols]),1)
roc_fpr_all['std'] = roc_fpr_all.apply(lambda x: np.std(x[cols]),1)

roc_fpr_all.to_csv('./'+args.prefix+'_results/roc_false_positive_rate.tsv', sep='\t')

cols = roc_tpr_all.columns
roc_tpr_all['mean'] = roc_tpr_all.apply(lambda x: np.mean(x[cols]),1)
roc_tpr_all['std'] = roc_tpr_all.apply(lambda x: np.std(x[cols]),1)

roc_tpr_all.to_csv('./'+args.prefix+'_results/roc_true_positive_rate.tsv', sep='\t')

cols = roc_auc_all.columns
roc_auc_all['mean'] = roc_auc_all.apply(lambda x: np.mean(x[cols]),1)
roc_auc_all['std'] = roc_auc_all.apply(lambda x: np.std(x[cols]),1)

roc_auc_all.to_csv('./'+args.prefix+'_results/roc_area_under_curve.tsv', sep='\t')

plt.figure(figsize=(6, 6))
plt.subplot(1, 1, 1)

fpr_err = np.array(roc_fpr_all['std'])
tpr_err = np.array(roc_tpr_all['std'])
m_auc = str(round(np.array(roc_auc_all['mean'])[0],3))
m_auc_std = str(round(np.array(roc_auc_all['std'])[0],3))

plt.plot(roc_fpr_all['mean'], roc_tpr_all['mean'], '-', color='black',
         label='Mean ROC curve\n(AUC='+m_auc+', STD='+m_auc_std+')')

#plt.fill_between(roc_fpr_all['mean'], 
#                 roc_tpr_all['mean']+tpr_err, roc_tpr_all['mean']-tpr_err,
#                 color='blue', alpha=0.2)

#plt.fill_betweenx(roc_tpr_all['mean'], 
#                  roc_fpr_all['mean']+fpr_err, roc_fpr_all['mean']-fpr_err,
#                  color='red', alpha=0.2)

plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate', fontsize=14)
plt.ylabel('True Positive Rate', fontsize=14)
plt.title('Mean receiver operating characteristic', fontsize=14)
plt.legend(loc="lower right")
plt.savefig('./'+args.prefix+'_results/classifier_plots/mean_roc.pdf',dpi=300, bbox_inches = 'tight')
plt.close()

print "Plotting mean Precision-Recall curve and error region"
print
cols = pre_all.columns
pre_all['mean'] = pre_all.apply(lambda x: np.mean(x[cols]),1)
pre_all['std'] = pre_all.apply(lambda x: np.std(x[cols]),1)

pre_all.to_csv('./'+args.prefix+'_results/precision.tsv', sep='\t')

cols = rec_all.columns
rec_all['mean'] = rec_all.apply(lambda x: np.mean(x[cols]),1)
rec_all['std'] = rec_all.apply(lambda x: np.std(x[cols]),1)

rec_all.to_csv('./'+args.prefix+'_results/recall.tsv', sep='\t')

cols = roc_auc_all.columns
av_pre_all['mean'] = av_pre_all.apply(lambda x: np.mean(x[cols]),1)
av_pre_all['std'] = av_pre_all.apply(lambda x: np.std(x[cols]),1)

av_pre_all.to_csv('./'+args.prefix+'_results/precision_recall_area_under_curve.tsv', sep='\t')

plt.figure(figsize=(6, 6))
plt.subplot(1, 1, 1)

pre_err = np.array(pre_all['std'])
rec_err = np.array(rec_all['std'])
m_auc = str(round(np.array(av_pre_all['mean'])[0],3))
m_auc_std = str(round(np.array(av_pre_all['std'])[0],3))

plt.plot(rec_all['mean'], pre_all['mean'], '-', color='black',
         label='Mean Precision-Recall curve\n(AUC='+m_auc+', STD='+m_auc_std+')')

#plt.fill_between(rec_all['mean'], 
#                 pre_all['mean']+pre_err, pre_all['mean']-pre_err,
#                 color='blue', alpha=0.2)

#plt.fill_betweenx(pre_all['mean'], 
#                  rec_all['mean']+rec_err, rec_all['mean']-rec_err,
#                  color='red', alpha=0.2)

plt.xlim([-0.05, 1.05])
plt.ylim(ymax=1.05)
plt.xlabel('Recall', fontsize=14)
plt.ylabel('Precision', fontsize=14)
plt.title('Precision-Recall', fontsize=14)
plt.legend(loc="lower left")
plt.savefig('./'+args.prefix+'_results/classifier_plots/mean_precision_recall.pdf',dpi=300, bbox_inches = 'tight')
plt.close()

print "Zipping folders from every run in a single file"
print 

Popen('tar -zcf ./'+args.prefix+'_results/classifier_plots/run_files.tar.gz ./'+args.prefix+'_results/classifier_plots/run_*/', shell=True)

print "Analysis complete. Thanks for using BioFeatureFinder."
print
sys.exit()
