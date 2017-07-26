## Load the packages required

import sys
import os
import argparse
from subprocess import PIPE, Popen
from multiprocessing.pool import ThreadPool
from multiprocessing.pool import Pool
import multiprocessing as mp
import warnings
import glob
from io import StringIO
import shutil
import tarfile
import pandas as pd
import numpy as np
import pybedtools
from pybedtools import BedTool
import pysam
import itertools

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
            "Check the help below and try to fix the arguments. If the error persists, please contact the corresponding author")
        print
        self.print_help()
        sys.exit(2)


## Assign input data as system variables

parser = MyParser(description='')

parser.add_argument('-b', '--bed', dest="bed_file",
                    help="List of BED intervals for creating the generic matrix.",
                    required=True)

parser.add_argument('-gen', '--genome', dest="genome_file",
                    help="Genome FASTA associated with the build of the GTF file",
                    required=True)

parser.add_argument('-o', '--outfile', dest="prefix",
                    help="prefix for use on the output files",
                    metavar="prefix", required=True)

parser.add_argument('-g', '--gtf', dest="gtf_file",
                    help="GTF file downloaded from Ensembl database. Available at http://www.ensembl.org/info/data/ftp/index.html. Used for analysis considering intronic, exonic, UTR and CDS.",
                    required=False)

parser.add_argument('-cs', '--conservation-scores', nargs='+', dest="con_files",
                    default=False,
                    help="bigWig file with phastCon scores for multiple alignments. Used as a measure of conservation of the features among the aligned species, obtained from 'http://hgdownload.cse.ucsc.edu/downloads.html' under 'Conservation scores' and downloading bigWig (.bw) files. Can take multiple files as input and accepts wildcard characters (*). REQUIRES bigWigAverageOverBed tool to be installed and available on PATH (can be obtained at UCSCs binaries directory (http://hgdownload.cse.ucsc.edu/admin/exe/). If no bigWig file is available, you can download the raw phastCon scores (.pp files) and create your own bigWig files using the wigToBigWig tool from the same repository. Default: False",
                    metavar="sp.phastCons*.bw", required=False)

parser.add_argument('-var', '--variation', nargs='+', dest="var_files", default=False,
                    help="Annotation file containing variation regions found in the genome (can be SNPs, strucutural variations, mutations or custom annotations). Can be obtained from UCSCs database or from Ensembl's GVF ftp directory. Can take multiple files as input and accepts wildcard characters (*). Default: False",
                    required=False)

parser.add_argument('-custom_seq', '--custom-sequence',
                    dest="custom_seq_ex", default=False,
                    help="Txt list containing custom sequences to be counted in the exons. Default: False",
                    metavar="custom_exon.txt", required=False)

parser.add_argument('-k', '--kmer', nargs='+', dest="kmer_list",
                    help="List of INT to create k-mers for counting. Default: False",
                    type=int, default=False, required=False)

parser.add_argument('-n', '--n_rand', dest="n_rand",
                    help="Number of times to shuffle BED input for background generator. Default: 3",
                    type=int, default=3, required=False)

parser.add_argument('-nuccont', '--nucleotide_content', dest="nuc_info",
                    default=1, metavar="nuc", required=False,
                    help="Defines the ammount of information included from the nucleotide sequence, 3 options available: 'Simple','Intermediate','Full'. Options:1 = Simple:[Length and pGC], 2 = Intermediate:[Length, pGC, pG, pC, pA, pT], 3 = Full:[All data from BedTools nucleotide sequence].' Default: 1 (Intermediate); p = percentage")

parser.add_argument("--keepBED", dest="keep_bed",
                    action="store_true", default=False,
                    help="Save the bed files generated for each class. Default: False")

parser.add_argument("--keepTEMP", dest="keep_temp",
                    action="store_true", default=False,
                    help="Keep the temporary files generated in each step. Default: False")

parser.add_argument("--ncores", dest="ncores", default=(mp.cpu_count() - 1),
                    help="Number of CPU cores used to process downtream/upstream exon search. Default:(ALL_CORES)-1",
                    type=int, metavar='INT')

parser.add_argument("--useThreads", dest="use_threads",
                    action="store_true", default=False,
                    help="Spawns multiple threads instead of multiple processes for parallelization. Is slower and uses mor RAM, but requires less CPU power. Default: False")

parser.add_argument("--debug", dest="debug", metavar='INT',
                    type=int, default=False,
                    help="Only seaches for the first N entries. Useful for debugin and checking if the code is working properly. Default: False")

args = parser.parse_args()


def concatFiles(input_files, output_file):
    files = glob.glob(input_files)
    with open(output_file, "wb") as fo:
        for f in sorted(files):
            with open(f, "rb") as fi:
                shutil.copyfileobj(fi, fo)

def save_bed(df, filename):
    df.to_csv(filename, compression='gzip', sep='\t', index=False,
              header=False)

def nuc_cont(bedtool):
    nuccont = bedtool.nucleotide_content(genome_fasta, s=True, seq=True)
    nuccont_df = pd.concat(nuccont.to_dataframe(
        names=['seqname', 'start', 'end', 'name', 'score', 'strand',
               '%AT', '%GC', '%A', '%C', '%G', '%T',
               '%N', '%O', 'length', 'seq'],
        iterator=True, chunksize=10000), ignore_index=True).ix[1:]
    nuccont_df['%AT'] = nuccont_df['%AT'].astype(float).round(5)
    nuccont_df['%GC'] = nuccont_df['%GC'].astype(float).round(5)
    nuccont_df['%A'] = nuccont_df.apply(
        lambda x: (float(x['%A']) / float(x['length'])), 1).round(5)
    nuccont_df['%C'] = nuccont_df.apply(
        lambda x: (float(x['%C']) / float(x['length'])), 1).round(5)
    nuccont_df['%G'] = nuccont_df.apply(
        lambda x: (float(x['%G']) / float(x['length'])), 1).round(5)
    nuccont_df['%T'] = nuccont_df.apply(
        lambda x: (float(x['%T']) / float(x['length'])), 1).round(5)
    nuccont_df['%N'] = nuccont_df.apply(
        lambda x: (float(x['%N']) / float(x['length'])), 1).round(5)
    nuccont_df['%O'] = nuccont_df.apply(
        lambda x: (float(x['%O']) / float(x['length'])), 1).round(5)

    if args.nuc_info == 1:
        nuccont_df = nuccont_df[['name', 'length', '%GC', 'seq']]
    elif args.nuc_info == 2:
        nuccont_df = nuccont_df[['name', 'length', '%GC', '%G', '%C', '%A', '%T', 'seq']]
    elif args.nuc_info == 3:
        pass

    return nuccont_df

def get_conservation_scores(con_file, df, cmd="bigWigAverageOverBed"):
    print
    print("Getting conservation from: " + str(con_file))
    print
    BedTool.from_dataframe(df).saveas('conservation_temp.bed')
    composed_command = " ".join([cmd, con_file, 'conservation_temp.bed',
                                 'conservation_result_temp.tab'])
    p = Popen(composed_command, shell=True)
    stdout, stderr = p.communicate()
    source = str(con_file.split('/')[-1])
    result = pd.concat(pd.read_table('conservation_result_temp.tab',
                                     names=('name', 'length', 'covered',
                                            'sum', 'mean_' + source,
                                            'mean0_' + source), iterator=True,
                                     chunksize=10000
                                     ), ignore_index=True)
    os.remove('conservation_temp.bed')
    os.remove('conservation_result_temp.tab')
    return result[['name', 'mean_' + source]]

def get_var_counts(bedtool, var_file):
    print
    print("Getting variation from: " + str(var_file))
    print
    var = BedTool(var_file)
    source = str(var_file).split('/')[-1]
    var_counts = pd.concat(
        bedtool.intersect(var, s=True, c=True).to_dataframe(iterator=True,
                                                            chunksize=10000
                                                            ),
        ignore_index=True).rename(
        columns={'thickStart': 'var_count_' + source})
    return var_counts[['name', 'var_count_' + source]]

def filter_columns(df):
    column_list = df.columns
    for i in range(len(column_list)):
        if df[column_list[i]].value_counts().shape[0] == 1:
            df.drop(column_list[i], 1, inplace=True)
        else:
            pass
    return df

def run_emboss_wordcount(filename, kmers):
    composed_command = " ".join(['wordcount -sequence emboss/fastas/'+filename, \
                                            '-wordsize '+str(kmers), \
                                            '-mincount 1',\
                                            '-outfile emboss/wordcount/'+filename+'.'+str(kmers)+'.wc',
                                            '2>>/dev/null'])
    p = Popen(composed_command, stdout=PIPE, shell=True)
    p.communicate()

def get_kmer_counts(df, df2, kmer_list):
    Popen('mkdir -p ./emboss', shell=True)
    Popen('mkdir -p ./emboss/fastas', shell=True)
    Popen('mkdir -p ./emboss/wordcount', shell=True)
    for i in range(len(df)):
        b = BedTool.from_dataframe(df.iloc[[i]])
        b.sequence(fi=genome_fasta, s=True, name=True,
                   fo='./emboss/fastas/'+df.iloc[[i]]['name'].to_string().split('range_id_')[1]+'.fa')

    filename = pd.DataFrame(glob.glob('./emboss/fastas/*'))[0].apply(lambda x: x.split('/')[-1])
    
    output = df2.copy()
    
    for i in range(len(kmer_list)):
        for j in range(len(filename)):
            run_emboss_wordcount(filename[j], 
                                 kmer_list[i])
       
    wc_list = glob.glob('./emboss/wordcount/*.wc')

    df_list = []

    for i in range(len(wc_list)):
        name = wc_list[i].split('/')[-1]
        name = name.split('.fa')[0]
        df = pd.read_table(wc_list[i],
                           index_col=0,
                           names=['range_id_'+name])
        df_list.append(df)

    kmer_df = pd.concat(df_list, axis=1, join='outer')
        
    kmer_df = kmer_df.T.add_prefix('kmer_count_').reset_index().fillna(0).rename(columns={'index':'name'})
    output = output.merge(kmer_df, on='name')
    return output

def get_chromsizes(genome_fasta, cmd="cut -f 1,2"):
    composed_command = " ".join([cmd, genome_fasta+'.fai', '>', genome_fasta+'.chromsizes'])
    p = Popen(composed_command, stdout=PIPE, shell=True)
    p.communicate()
    
def shuffle_bedtools_gtf(number):
    composed_command = " ".join(['bedtools shuffle', '-i', args.bed_file, 
                                 '-g', args.genome_file+'.chromsizes', 
                                 '-incl', args.gtf_file,
                                 '-chromFirst',
                                 '> shuffled.entry.'+str(number)+'.bed'])
    p = Popen(composed_command, stdout=PIPE, shell=True)
    p.communicate()    
    
def shuffle_bedtools_no_gtf(number):
    composed_command = " ".join(['bedtools shuffle', '-i', args.bed_file, 
                                 '-g', args.genome_fasta+'.chromsizes',
                                 '-chromFirst',
                                 '> shuffled.entry.'+str(number)+'.bed'])
    p = Popen(composed_command, stdout=PIPE, shell=True)
    p.communicate()

def get_data(df, name, matrix):
    print
    print("Starting " + name)

    bedtool = BedTool.from_dataframe(df).saveas()
    a = nuc_cont(bedtool)
    
    if args.kmer_list:
        print
        print ("Starting K-mer counting")
        a = get_kmer_counts(df, a, args.kmer_list)
    elif not args.kmer_list:
        pass

    if args.var_files:
        var_files_list = glob.glob(str(args.var_files))
        for i in range(len(var_files_list)):
            var = get_var_counts(bedtool, var_files_list[i])
            a = a.merge(var, on='name')
    elif not args.var_files:
        pass

    if args.con_files:
        con_files_list = glob.glob(str(args.con_files))
        for i in range(len(con_files_list)):
            con = get_conservation_scores(con_files_list[i], df)
            a = a.merge(con, on='name')
    elif not args.con_files:
        pass

    if str(args.nuc_info) == 3:
        z = a.drop(['seqname', 'start', 'end', 'score', 'strand', 'seq'], 1)
    else:
        z = a.drop('seq', 1)

    z = z.set_index('name').add_suffix('_' + name).reset_index()
    z = filter_columns(z)
    print
    print(name + " finished")
    return matrix.merge(z, on='name').drop_duplicates()    

print
print "Starting datamatrix assembly process"

##Load the genome file that matches the version of the GTF you are using. Pysam will be used to build an index of
##the FASTA file.

print
print "Loading genome and creating FASTA index file and extracting chrom sizes"

genome_fasta = args.genome_file
pysam.faidx(genome_fasta)
get_chromsizes(args.genome_file)

## Load the reference GTF file and convert it to a DataFrame

if not args.debug:
    pass
else:
    print
    print("Running in debug mode. Only the first " + str(args.debug) + " entries will be used.")

if not args.use_threads:
    pass
else:
    print
    print("Running with multiple threads. Watch out for that RAM!")

if args.gtf_file:
    print
    print "Loading reference GTF file"
    gtf_ref = BedTool(args.gtf)   

print
print "Sorting input bed file."

input_bed = BedTool(args.bed_file).sort().saveas('input_list.bed')
    
if args.n_rand >= 1:
    print
    print "Shuffling input bed in the genome and generating randomized background."
    print
    print "Generating "+str(args.n_rand)+" times the size of input list for random background."
    n_rand = np.arange(args.n_rand)
        
    if args.gtf_file:
        if __name__ == '__main__':
            p = Pool((mp.cpu_count() - 1))
            p.map(shuffle_bedtools_gtf, n_rand)
       
    elif not args.gtf_file:
        if __name__ == '__main__':
            p = Pool((mp.cpu_count() - 1))
            p.map(shuffle_bedtools_no_gtf, n_rand)
        
    concatFiles('shuffled.entry.*.bed','shuffled.bed')
    list(map(os.remove, glob.glob("shuffled.entry.*.bed")))
        
    cat_command = " ".join(['cat shuffled.bed input_list.bed > genomic_ranges.bed'])
    p = Popen(cat_command, stdout=PIPE, shell=True)
    p.communicate()
    
if args.n_rand == 0:
    print
    print "Skipping randomized background step"
    cat_command = " ".join(['cat input_list.bed > genomic_ranges.bed'])
    p = Popen(cat_command, stdout=PIPE, shell=True)
    p.communicate()
        
    
bed = BedTool('genomic_ranges.bed').to_dataframe()
bed['name'] = 'range_id_R' + (bed.index + 1).astype(str) + '_' + \
                                bed['chrom'].astype(str) + '_' + \
                                bed['start'].astype(str) + '_' + \
                                bed['end'].astype(str) + '_' + \
                                bed['strand'].astype(str)
    
bed = bed[['chrom', 'start', 'end', 'name', 'score', 'strand']
         ].drop_duplicates().sort_values(by=['chrom', 'start'])
    
if not args.debug:
    pass
else:
    bed = bed.head(args.debug)
    
matrix_bed = pd.DataFrame(bed['name'])
matrix_bed = get_data(bed, 'genomic.ranges', matrix_bed)
    
matrix_bed.set_index('name').drop_duplicates().to_csv(str(args.prefix) + '.genomic.ranges.datamatrix.tsv', sep='\t')

bed.to_csv('./genomic_ranges.bed', header=False, index=False, sep='\t')

print
print("Data matrices build complete")
