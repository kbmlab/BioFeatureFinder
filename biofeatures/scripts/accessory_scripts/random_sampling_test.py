import sys
import time
import argparse
from subprocess import Popen
from itertools import cycle
from multiprocessing.pool import Pool
import multiprocessing as mp
import warnings
import glob

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
pd.options.mode.chained_assignment = None  # default='warn'

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        print()
        print("The following error ocurred in argument parsing:")
        sys.stderr.write('error: %s\n' % message)
        print()
        print("Check the help below and try to fix the arguments. If the error persists, please contact the corresponding author")
        print()
        self.print_help()
        sys.exit(2)

## Assign input data as system variables

parser = MyParser(description='')

parser.add_argument('-df', '--data_frame', dest="data", 
                  help="Data frame with features created by 'build_datamatrix.py'.", required=True)

parser.add_argument('-reps', '--repetitions', dest="nreps", default=100,
                  help="Number of random samplings to be performed in each feature. Default:100", required=False)

parser.add_argument('-filter', '--filter_columns', dest='filter_list', default=False, required=False,
                    help="Input table with strings to filter out columns. Can be used with partial string match.")

parser.add_argument("--ncores", dest="ncores", default=(mp.cpu_count()-1), help="Number of CPU cores used to process downtream/upstream exon search. Default:(ALL_CORES)-1", type=int, metavar='INT')

parser.add_argument('-o', '--outfile', dest="prefix", 
                  help="prefix for use on the output files", metavar="prefix", required=True)

args = parser.parse_args()

t0 = time.time()

##Load the data matrix in an object and set the name column as index

print()
print("Load data matrix")

matrix = pd.concat(pd.read_table(str(args.data),iterator=True, chunksize=10000), ignore_index=True).set_index('name')

##Filter or drop columns and set the name as index for the analysis
   
if not args.filter_list:
    pass
else:
    filter_list = list(pd.read_table(args.filter_list, header=None)[0])
    print()
    print("Filtering matrix for removing columns: "+str(filter_list))
    for i in range(len(filter_list)):
        matrix = matrix.drop(
            matrix.filter(like=str(filter_list[i]), axis=1).columns,1)
        
##Set the number of repetitions to the used and the fractions to be sampled

print()
print("Setting testing parameters")

reps = np.arange(int(args.nreps))
f_samp = list((0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

##Create the folder for saving the output files and graphs

print()
print("Creating output folders")

Popen('mkdir random_sampling_test', shell=True)
Popen('mkdir random_sampling_test/'+str(args.prefix), shell=True)
        
##Create a function for reading the features and applying statistical tests. The ouput will be CSV files with the data and graphs for the individual features.

def get_rand_sampling_stats(features):
    name = features.split('/')[-1]
    test_reps = pd.DataFrame()
    test_reps[str(name)] = reps+1
    for j in range(len(f_samp)):
        test_reps['ks_'+str(int(f_samp[j]*100))] = test_reps[str(name)
                                                            ].apply(lambda x: stats.ks_2samp(
                matrix[str(features)].astype(float),
                matrix.sample(frac=f_samp[j])[str(features)].astype(float))[1])
        
    Tr = test_reps.set_index(str(name)).T
    Tr['median'] = Tr.apply(lambda x: np.median(x),1)
    Tr['std'] = Tr.apply(lambda x: np.std(x),1)
    Tr.T.to_csv('random_sampling_test/'+str(args.prefix)+'/'+str(name)+'_test_for_random_sampling.csv')
    plt.errorbar(y=Tr['median'], x=f_samp, yerr=Tr['std'], fmt=':',)
    plt.title(str(name))
    plt.xlabel('fraction of exons sampled')
    plt.ylabel('ks p-val')
    plt.xlim(0,1.1)
    plt.savefig('random_sampling_test/'+str(args.prefix)+'/'+str(name)+'_test_for_random_sampling.pdf',
               dpi=300,bbox_inches='tight')

features = matrix.columns    

print()
print("Running tests with "+str(args.nreps)+" repetitions per feature and saving output files")

if __name__ == '__main__':
    p = Pool(args.ncores)
    p.map(get_rand_sampling_stats, features)

print()
print("Plotting composed figures for all features")
    
files = glob.glob('random_sampling_test/'+str(args.prefix)+'/*.csv')
cycol = cycle('bgrcmk').__next__

for i in range(len(files)):
    df = pd.read_csv(files[i], index_col=0).T
    
    coefficients = np.polyfit(np.log(f_samp),df['median'],2) # Use log(x) as the input to polyfit.
    fit = np.poly1d(coefficients) 
    
    color = cycol()
    plt.errorbar(y=df['median'], x=f_samp, yerr=df['std'], fmt=':', c=color)
    plt.plot(f_samp,fit(np.log(f_samp)),"--", label="fit", c=color)
    plt.xlabel('fraction of total exons sampled')
    plt.ylabel('median ks p-val')
    plt.xlim(0,1.1)
#    plt.ylim(ymin=0)

plt.axhline(0.95)
plt.xlim(0,1.1)
#plt.ylim(ymin=0)
plt.title('Random sampling fractions for '+str(args.prefix)+' exons features')
plt.savefig('random_sampling_test/'+str(args.prefix)+'_random_sampling.pdf',
            dpi=300,bbox_inches='tight')

t1 = time.time()

print() 
print("Done! Running time: "+str(t1-t0)+" seconds.")