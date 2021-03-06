#!/usr/bin/env python
## Load the packages required

import sys
import os
import argparse
from subprocess import PIPE, Popen
from multiprocessing.pool import Pool
import multiprocessing as mp
import warnings
import glob
import shutil
import pandas as pd
import numpy as np
import pybedtools
from pybedtools import BedTool
import pysam
import itertools

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
pd.options.mode.chained_assignment = None  # default='warn'

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
        print()
        self.print_help()
        print()
        print(gandalf)
        print()    
        print()
        print("The following error ocurred in argument parsing:")
        sys.stderr.write('error: %s\n' % message)
        print()
        print(
            "Check the help and try to fix the arguments. If the error persists, please contact the corresponding author")
        print()   
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

parser.add_argument('-i', '--input', dest="input_file", 
                    help="Input list of intervals to analyze.",
                    required=True)

parser.add_argument('-gen', '--genome', dest="genome_file",
                    help="Genome FASTA associated with the interval file.",
                    required=True)

parser.add_argument('-o', '--outfile', dest="outfile", type=str,
                    help="Name for output files generated.",
                    required=True)

parser.add_argument('-g', '--gtf', dest="gtf_file",
                    help="GTF file for guided background generation. Default: False",
                    required=False)
                    
parser.add_argument('-nuc', '--nucleotide_content', dest="nuc_info",
                    default='simple', metavar="nuc", required=False, type=str,
                    help = "Defines the ammount of information included from the nucleotide sequence, 3 options available: simple [Length and GCp], intermediate [Length, GCp, Gp, Cp, Ap, Tp], full [All data from BedTools nucleotide sequence]. Default: simple; p = percentage")
 
parser.add_argument('-bw', '--bigwig-scores', nargs='+', dest="bw_files",
                    default=False,
                    help="This option takes as input bigWig files with scores and calculates the mean score over each input region. REQUIRES bigWigAverageOverBed tool to be installed and available on PATH (can be obtained at UCSCs binaries directory (http://hgdownload.cse.ucsc.edu/admin/exe/).",
                    metavar="somefile.bw", required=False)

parser.add_argument('-var', '--variation', nargs='+', dest="var_files", default=False,
                    help="This option takes as input BED files containing regions of biological featues found in the genome (can be SNPs, strucutural variations, mutations or custom annotations) and counts the number of occurences in each input region. Can take multiple files as input and accepts wildcard characters (*). Default: False",
                    metavar='somefile.bed', required=False)

parser.add_argument('-f','--fasta', dest="create_fastas", action='store_true',
                    help="Use this option to create fasta files for each entry in the analysis. Required for k-mer search (-k) and structural MFE (-s and -sg). Default: False",
                    required=False, default=False)

parser.add_argument('-k', '--kmer', nargs='+', dest="kmer_list",
                    help="List of INT to create k-mers for counting. Default: False",
                    type=int, default=False, required=False)

parser.add_argument('-rf','--rnafold', dest="rnafold",
                    action="store_true", default=False, required=False,
                    help="Run RNAFold (from Vienna RNA package) on each entry and extract MFE values. Requires Vienna RNA Package installed locally (https://www.tbi.univie.ac.at/RNA/) and available on PATH. Default: False")

parser.add_argument('-qg','--qgrs', dest="qgrs_mapper",
                    action="store_true", default=False, required=False,
                    help="Run QGRS Mapper on each entry and extract G-Quadruplex scores. Requires QGRS Mapper installed locally (https://github.com/freezer333/qgrs-cpp; http://bioinformatics.ramapo.edu/temp/qgrs/index.php) and available on PATH. Default: False")

parser.add_argument('-custom_seq', '--custom-sequence',
                    dest="custom_seq_ex", default=False,
                    help="Txt list containing custom sequences to be counted in the exons. Default: False",
                    metavar="custom_exon.txt", required=False)

parser.add_argument('-n', '--n_rand', dest="n_rand",
                    help="Number of times to shuffle BED input for background generator. Default: 3",
                    type=int, default=3, required=False)

parser.add_argument("--keepBED", dest="keep_bed",
                    action="store_true", default=False,
                    help="Save the bed files generated for each class. Default: False")

parser.add_argument("--keepTEMP", dest="keep_temp",
                    action="store_true", default=False,
                    help="Keep the temporary files generated in each step. Default: False")

parser.add_argument('-c',"--ncores", dest="ncores", default=(mp.cpu_count() - 1),
                    help="Number of CPU cores used to process downtream/upstream exon search. Default:(ALL_CORES)-1",
                    type=int, metavar='INT')

parser.add_argument("-d","--debug", dest="debug", metavar='INT',
                    type=int, default=False,
                    help="Only seaches for the first N entries. Useful for debugin and checking if the code is working properly. Default: False")

args = parser.parse_args()

if args.kmer_list is True and args.create_fastas is False:
    parser.error("--kmer requires --fasta.")
    
if args.rnafold is True and args.create_fastas is False:
    parser.error("--rnafold requires --fasta.")
       
def nuc_cont(bedtool):
    nuccont = bedtool.nucleotide_content(args.genome_file, s=strd, seq=True)
    
    if strd == False:
        nuccont_df = pd.concat(nuccont.to_dataframe(
                names=['seqname', 'start', 'end', 'name', 'score',
                       '%AT', '%GC', '%A', '%C', '%G', '%T',
                       '%N', '%O', 'length', 'seq'],
                iterator=True, chunksize=10000), 
        ignore_index=True, sort=False).ix[1:]
    
    elif strd == True:
        nuccont_df = pd.concat(nuccont.to_dataframe(
            names=['seqname', 'start', 'end', 'name', 'score', 'strand',
                   '%AT', '%GC', '%A', '%C', '%G', '%T',
                   '%N', '%O', 'length', 'seq'],
            iterator=True, chunksize=10000), 
        ignore_index=True, sort=False).ix[1:]
    
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

    if args.nuc_info == 'simple':
        nuccont_df = nuccont_df[['name', '%GC', 'seq']]
    elif args.nuc_info == 'intermediate':
        nuccont_df = nuccont_df[['name', '%GC', '%G', '%C', '%A', '%T', 'seq']]
    elif args.nuc_info == 'full':
        pass

    return nuccont_df

def get_bigwig_scores(bw_file, df, cmd="bigWigAverageOverBed"):
    print()
    print(("Getting average scores from: " + str(bw_file)))
    print()
    BedTool.from_dataframe(df).saveas(args.outfile+'.datamatrix/bigwig_temp.bed')
    composed_command = " ".join([cmd, bw_file, args.outfile+'.datamatrix/bigwig_temp.bed',
                                 args.outfile+'.datamatrix/bigwig_result_temp.tab'])
    p = Popen(composed_command, shell=True)
    stdout, stderr = p.communicate()
    source = str(bw_file.split('/')[-1])
    result = pd.concat(pd.read_table(args.outfile+'.datamatrix/bigwig_result_temp.tab',
                                     names=('name', 'length', 'covered',
                                            'sum', 'mean_' + source,
                                            'mean0_' + source), iterator=True,
                                     chunksize=10000
                                     ), ignore_index=True, sort=False)
    os.remove(args.outfile+'.datamatrix/bigwig_temp.bed')
    os.remove(args.outfile+'.datamatrix/bigwig_result_temp.tab')
    return result[['name', 'mean_' + source]]

def get_var_counts(bedtool, var_file):
    print()
    print(("Getting variation from: " + str(var_file)))
    print()
    var = BedTool(var_file).sort().remove_invalid().saveas(args.outfile+'.datamatrix/varfile')#.sort()
    source = str(var_file).split('/')[-1]
    
    if strd == True:
        
        if 'strand' not in list(BedTool(var[0:1]).saveas().to_dataframe().columns):
            
            print(source+' does not contain +/- strand information. Running in unstranded mode')
            print()
            var_counts = pd.concat(
                    bedtool.intersect(var, 
                                      s=False, 
                                      c=True, 
                                      sorted=False
                                      ).to_dataframe(iterator=True,
                                                     chunksize=10000),
                                   ignore_index=True, 
                                   sort=False)    
        else:
            var_counts = pd.concat(
                    bedtool.intersect(var, 
                                      s=strd, 
                                      c=True, 
                                      sorted=False
                                      ).to_dataframe(iterator=True,
                                                     chunksize=10000),
                                   ignore_index=True, 
                                   sort=False)
    
    elif strd == False:
        var_counts = pd.concat(
        bedtool.intersect(var, 
                          s=False, 
                          c=True, 
                          sorted=False
                          ).to_dataframe(iterator=True,
                                         chunksize=10000),
                       ignore_index=True, 
                       sort=False)  
        
    var_counts = var_counts.rename(
            columns={var_counts.columns[-1]: 'var_count_' + source})    
    return var_counts[['name', 'var_count_' + source]]

def filter_columns(df):
    column_list = df.columns
    for i in range(len(column_list)):
        if df[column_list[i]].value_counts().shape[0] == 1:
            df.drop(column_list[i], 1, inplace=True)
        else:
            pass
    return df

def make_fasta(entry):
    df_slice = bed.iloc[[entry]]
    b = BedTool.from_dataframe(df_slice)
    b.sequence(fi=args.genome_file, s=strd, name=True,
               fo=args.outfile+'.datamatrix/temp/fastas/'+df_slice['name'].to_string().split('range_id_')[1]+'.fa')
    
def run_emboss_wordcount(arg):
    file_entry = arg[0]
    kmer = arg[1]
    composed_command = " ".join(['wordcount -sequence '+args.outfile+'.datamatrix/temp/fastas/'+file_entry, \
                                            '-wordsize '+str(kmer), \
                                            '-mincount 1',\
                                            '-outfile '+args.outfile+'.datamatrix/temp/emboss/'+file_entry+'.'+str(kmer)+'.wc',
                                            '2>>/dev/null'])
    p = Popen(composed_command, stdout=PIPE, shell=True)
    p.communicate()

def get_kmer_counts(kmer_list):
    Popen('mkdir -p '+args.outfile+'.datamatrix/temp/emboss', shell=True)
    filename = pd.DataFrame(glob.glob(args.outfile+'.datamatrix/temp/fastas/*'))[0].apply(lambda x: x.split('/')[-1])

    for i in range(len(args.kmer_list)):
        
        kmer = args.kmer_list[i]
        
        if __name__ == '__main__':
            p = Pool(args.ncores)
            p.map(run_emboss_wordcount, zip(filename, itertools.repeat(args.kmer_list[i])))
       
        wc_list = glob.glob(args.outfile+'.datamatrix/temp/emboss/*.'+str(kmer)+'.wc')
        df_list = []

        for i in range(len(wc_list)):
            name = wc_list[i].split('/')[-1]
            name = name.split('.fa')[0]
            df = pd.read_table(wc_list[i],
                               index_col=0,
                               names=['range_id_'+name])
            df_list.append(df)

        #with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        #    kmer_df = executor.submit(pd.concat, df_list, axis=1, join='outer').result()
        
        kmer_df = pd.concat(df_list, axis=1, join='outer', sort=False)
        
        kmer_df = kmer_df.T.add_prefix('kmer_count_').reset_index()
        kmer_df = kmer_df.groupby('index').sum()
        kmer_df = kmer_df.reset_index().fillna(0).rename(columns={'index':'name'})
        kmer_df.to_csv(args.outfile+".datamatrix/temp/data_kmer_"+str(kmer)+"_results.csv", index=False)
        del kmer_df
    #df2 = df2.merge(kmer_df, on='name')
    #return df2

def run_rnafold(files):
    p = Popen("RNAfold -i "+args.outfile+".datamatrix/temp/fastas/"+files+\
              " --gquad --noPS | sed -n 3p | rev | cut -f 1 -d ' ' | rev | sed 's/(//g'| sed 's/)//g' > "+args.outfile+\
              ".datamatrix/temp/rnafold/"+files+\
              ".MFE", shell=True)
    p.communicate()

def get_MFE_scores():
    Popen('mkdir -p '+args.outfile+'.datamatrix/temp/rnafold', shell=True)
    fasta_files = pd.DataFrame(glob.glob(args.outfile+'.datamatrix/temp/fastas/*'))[0].apply(lambda x: x.split('/')[-1])
       
    if __name__ == '__main__':
        p = Pool(args.ncores)
        p.map(run_rnafold, fasta_files)
        
    rf_list = glob.glob(args.outfile+'.datamatrix/temp/rnafold/*.MFE')

    df_list = []

    for i in range(len(rf_list)):
        name = rf_list[i].split('/')[-1]
        name = name.split('.fa')[0]
        df = pd.read_table(rf_list[i],
                           names=['range_id_'+name])
        df_list.append(df)
    
    #with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
    #    mfe_df = executor.submit(pd.concat, df_list, axis=1, join='outer').result()
    
    mfe_df = pd.concat(df_list, axis=1, join='outer', sort=False)
    
    mfe_df = mfe_df.T.rename(columns={0:'MFE'}).reset_index().fillna(0).rename(columns={'index':'name'})
    mfe_df.to_csv(args.outfile+".datamatrix/temp/data_rnafold_results.csv", index=False)
    del mfe_df
    #df_input = df_input.merge(mfe_df, on='name')
    #return df_input

def run_qgrs_mapper(files):
    p = Popen('qgrs -csv -i '+args.outfile+'.datamatrix/temp/fastas/'+files+' -o '+args.outfile+'.datamatrix/temp/qgrs/'+files+'.qgrs.csv 2>>/dev/null', shell=True)
    p.communicate()
    
def get_QGRS_scores():
    Popen('mkdir -p '+args.outfile+'.datamatrix/temp/qgrs', shell=True)
    fasta_files = pd.DataFrame(glob.glob(args.outfile+'.datamatrix/temp/fastas/*'))[0].apply(lambda x: x.split('/')[-1])
       
    if __name__ == '__main__':
        p = Pool(args.ncores)
        p.map(run_qgrs_mapper, fasta_files)
        
    qgrs_list = glob.glob(args.outfile+'.datamatrix/temp/qgrs/*.qgrs.csv')

    df_list = []

    for i in range(len(qgrs_list)):
        name = qgrs_list[i].split('/')[-1]
        name = name.split('.fa.qgrs')[0]
        df = pd.read_csv(qgrs_list[i], skiprows=1)
        df = df.rename(columns={'GS':'range_id_'+name})
        df = df[['range_id_'+name]]
        df_list.append(df)

    #with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
    #    qgrs_df = executor.submit(pd.concat, df_list, axis=1, join='outer').result()    
    
    qgrs_df = pd.DataFrame(pd.concat(df_list, axis=1, join='outer', sort=False).max())
    
    qgrs_df = qgrs_df.rename(columns={0:'max_QGRS_score'}).reset_index().fillna(0).rename(columns={'index':'name'})
    qgrs_df.to_csv(args.outfile+".datamatrix/temp/data_qgrs_results.csv", index=False)
    del qgrs_df
    #df_input = df_input.merge(qgrs_df, on='name', how='outer')     
    #return df_input

def get_chromsizes(genome_file, cmd="cut -f 1,2"):
    composed_command = " ".join([cmd, genome_file+'.fai', '>', genome_file+'.chromsizes'])
    p = Popen(composed_command, stdout=PIPE, shell=True)
    p.communicate()
    
def shuffle_bedtools_gtf(number):
    composed_command = " ".join(['bedtools shuffle', '-i', args.input_file, 
                                 '-g', args.genome_file+'.chromsizes', 
                                 '-incl', args.gtf_file,
                                 '-chromFirst',
                                 '> '+args.outfile+'.datamatrix/shuffled.entry.'+str(number)+'.bed'])
    p = Popen(composed_command, stdout=PIPE, shell=True)
    p.communicate()    
    
def shuffle_bedtools_no_gtf(number):
    composed_command = " ".join(['bedtools shuffle', '-i', args.input_file, 
                                 '-g', args.genome_file+'.chromsizes',
                                 '-chromFirst',
                                 '> '+args.outfile+'.datamatrix/shuffled.entry.'+str(number)+'.bed'])
    p = Popen(composed_command, stdout=PIPE, shell=True)
    p.communicate()
    
def concatFiles(input_files, output_file):
    files = glob.glob(input_files)
    with open(output_file, "wb") as fo:
        for f in sorted(files):
            with open(f, "rb") as fi:
                shutil.copyfileobj(fi, fo)

def save_bed(df, filename):
    df.to_csv(filename, compression='gzip', sep='\t', index=False,
              header=False)

def get_data(df):

    Popen('mkdir -p ./'+args.outfile+".datamatrix/temp/", shell=True)
    
    bedtool = BedTool.from_dataframe(df).sort().saveas(args.outfile+'.datamatrix/bedtool_df.bed')
    a = nuc_cont(bedtool)   
    
    if args.var_files:
        var_files_list = list(args.var_files)
        for i in range(len(var_files_list)):
            var = get_var_counts(bedtool, var_files_list[i])
            a = a.merge(var, on='name')
    elif not args.var_files:
        pass

    if args.bw_files:
        bw_files_list = list(args.bw_files)
        for i in range(len(bw_files_list)):
            bw = get_bigwig_scores(bw_files_list[i], df)
            a = a.merge(bw, on='name')
    elif not args.con_files:
        pass

    if args.kmer_list:
        print()
        print("Starting K-mer counting")
        get_kmer_counts(args.kmer_list)
    elif not args.kmer_list:
        pass

    if args.rnafold == True:
        print()
        print("Starting RNAfold for MFE scoring")
        get_MFE_scores()
    elif args.rnafold == False:
        pass
    
    if args.qgrs_mapper == True:
        print()
        print("Starting QGRS Mapper for G-Quadruplex scoring")
        get_QGRS_scores()
    elif args.qgrs_mapper == False:
        pass        
 
    if str(args.nuc_info) == 'full':
        z = a.drop(['length','seqname', 'start', 'end', 
                    'score', 'strand', 'seq'], 1)
    else:
        z = a.drop('seq', 1)
    
    if (args.rnafold == True) or (args.qgrs_mapper == True) or (args.kmer_list is not False):
        z.to_csv(args.outfile+".datamatrix/temp/data_generic_results.csv", index=False)
        
        temp_files = glob.glob(args.outfile+".datamatrix/temp/*.csv")
        
        z_list = []
        
        for i in range(len(temp_files)):
            df = pd.read_csv(temp_files[i], index_col=0);
            if ((args.nuc_info == 'full') & ("data_generic_results.csv" in temp_files[i])):
                df = df.set_index('name')
            else:
                pass
            z_list.append(df)
        
        #with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        #    z = executor.submit(pd.concat, z_list, axis=1, join='outer').result()
        
        z = pd.concat(z_list, axis=1, join='outer', sort=False)
        z = z.reset_index().rename(columns={'index':'name'}).fillna(0)
    else:
        pass

    z = filter_columns(z)
    
    return z.drop_duplicates()    

if not args.debug:
    pass
else:
    print()
    print(("Running in debug mode. Only the first " + str(args.debug) + " entries will be used."))

print()
print("Starting datamatrix assembly process")

Popen('mkdir '+args.outfile+'.datamatrix', shell=True)

print()
print("Sorting input bed file.")

input_bed = BedTool(args.input_file).sort().saveas(args.outfile+'.datamatrix/input_list.bed')

if 'strand' in list(BedTool(input_bed[0:1]).saveas().to_dataframe().columns):
    print("Strand information found in input file. Running in stranded mode.")
    print()
    strd = True
else:
    print("Strand information NOT found in input file. Running in unstranded mode.")
    print()
    strd = False

##Load the genome file that matches the version of the GTF you are using. Pysam will be used to build an index of
##the FASTA file.

print()
print("Loading genome and checking for fasta index and chromsizes")

if os.path.isfile(args.genome_file+'.fai') == True:
    print()
    print("Genome FASTA index found.")
    pass
else:
    print()
    print("Genome FASTA index not found. Generating now...")
    pysam.faidx(args.genome_file)
    
if os.path.isfile(args.genome_file+'.chromsizes') == True:
    print()
    print("Genome chromsizes found.")
    pass
else:
    print()
    print("Genome chromsizes not found. Generating now...")
    get_chromsizes(args.genome_file)

## Load the reference GTF file and convert it to a DataFrame

if args.gtf_file:
    print()
    print("Loading reference GTF file")
    gtf_ref = BedTool(args.gtf_file)   
    
if args.n_rand >= 1:
    print()
    print("Shuffling input bed in the genome and generating randomized background.")
    print()
    print("Generating "+str(args.n_rand)+" times the size of input list for random background.")
    n_rand = np.arange(args.n_rand)
        
    if args.gtf_file:
        if __name__ == '__main__':
            p = Pool((args.ncores))
            p.map(shuffle_bedtools_gtf, n_rand)
       
    elif not args.gtf_file:
        if __name__ == '__main__':
            p = Pool((args.ncores))
            p.map(shuffle_bedtools_no_gtf, n_rand)
        
    concatFiles(args.outfile+'.datamatrix/shuffled.entry.*.bed',args.outfile+'.datamatrix/shuffled.bed')
    list(map(os.remove, glob.glob(args.outfile+".datamatrix/shuffled.entry.*.bed")))
        
    cat_command = " ".join(['cat '+args.outfile+'.datamatrix/shuffled.bed '+\
                            args.outfile+'.datamatrix/input_list.bed > '+\
                            args.outfile+'.datamatrix/genomic_ranges.bed'])
    p = Popen(cat_command, stdout=PIPE, shell=True)
    p.communicate()
    
if args.n_rand == 0:
    print()
    print("Skipping randomized background step")
    cat_command = " ".join(['cat '+args.outfile+'.datamatrix/input_list.bed > '+args.outfile+'.datamatrix/genomic_ranges.bed'])
    p = Popen(cat_command, stdout=PIPE, shell=True)
    p.communicate()
        
    
bed = BedTool(args.outfile+'.datamatrix/genomic_ranges.bed').to_dataframe().reset_index()

if strd == False:
    bed['name'] = 'range_id_R' + (bed['index'] + 1).astype(str) + '_' + \
                              bed['chrom'].astype(str) + '_' + \
                              bed['start'].astype(str) + '_' + \
                              bed['end'].astype(str)
    
    bed = bed[['chrom', 'start', 'end', 'name', 'score']
             ].drop_duplicates().sort_values(by=['chrom', 'start'])

elif strd == True:
    bed['name'] = 'range_id_R' + (bed['index'] + 1).astype(str) + '_' + \
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

BedTool.from_dataframe(bed).saveas(args.outfile+'.datamatrix/genomic_ranges.bed')

print()
print("Generating FASTA sequences for each entry")

if args.create_fastas == True:
    Popen('mkdir -p '+args.outfile+'.datamatrix/temp/fastas', shell=True)
    entry_list = list(range(len(bed)))

    if __name__ == '__main__':
        p = Pool((args.ncores))
        p.map(make_fasta, entry_list)
else:
    pass

print()
print("Starting data extraction and matrix assembly")

matrix_bed = get_data(bed)

print()
print("Saving final datamatrix in: "+args.outfile+'.datamatrix/'+str(args.outfile)+'.datamatrix.tsv')

matrix_bed.set_index('name').drop_duplicates().to_csv(
        args.outfile+'.datamatrix/'+str(args.outfile) + '.datamatrix.tsv', sep='\t')

print()
print("Cleaning up pybedtools temporary files")

pybedtools.helpers.cleanup(verbose=False, remove_all=False)

if args.keep_bed == True:
    Popen('rm '+args.outfile+'.datamatrix/varfile', 
          shell=True)
else:
    Popen('rm '+args.outfile+'.datamatrix/shuffled.bed '+\
                args.outfile+'.datamatrix/input_list.bed '+\
                args.outfile+'.datamatrix/varfile '+\
		args.outfile+'.datamatrix/bedtool_df.bed '+\
		args.outfile+'.datamatrix/genomic_ranges.bed',
          shell=True)
    
if args.keep_temp == True:
    pass
else:
    Popen('rm -r '+args.outfile+'.datamatrix/temp/', 
          shell=True)

print()
print("Thank you for using BioFeatureFinder")
print()
