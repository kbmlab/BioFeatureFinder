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

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

pd.options.mode.chained_assignment = None  # default='warn'


##Load the parser for arguments

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        print()
        print("The following error ocurred in argument parsing:")
        sys.stderr.write('error: %s\n' % message)
        print()
        print(
            "Check the help below and try to fix the arguments. If the error persists, please contact the corresponding author")
        print()
        self.print_help()
        sys.exit(2)


## Assign input data as system variables

parser = MyParser(description='')

parser.add_argument('-g', '--gtf', dest="gtf_file",
                    help="GTF file downloaded from Ensembl database. Available at http://www.ensembl.org/info/data/ftp/index.html",
                    required=True)

parser.add_argument('-gen', '--genome', dest="genome_file",
                    help="Genome FASTA associated with the build of the GTF file",
                    required=True)

parser.add_argument('-o', '--outfile', dest="prefix",
                    help="prefix for use on the output files",
                    metavar="prefix", required=True)

parser.add_argument('-cs', '--conservation_scores', dest="con_files",
                    default=False,
                    help="bigWig file with phastCon scores for multiple alignments. Used as a measure of conservation of the features among the aligned species, obtained from 'http://hgdownload.cse.ucsc.edu/downloads.html' under 'Conservation scores' and downloading bigWig (.bw) files. Can take multiple files as input by using wildcard characters (*) with single quotation marks (ex. 'hg38.*.bw'). REQUIRES bigWigAverageOverBed tool to be installed and available on PATH (can be obtained at UCSCs binaries directory (http://hgdownload.cse.ucsc.edu/admin/exe/). If no bigWig file is available, you can download the raw phastCon scores (.pp files) and create your own bigWig files using the wigToBigWig tool from the same repository. Default: False",
                    metavar="sp.phastCons*.bw", required=False)

parser.add_argument('-cpg', '--cpg_islands', dest="cpg_file", default=False,
                    help='BED file containing the CpG islands sites found in the genome. Can be obtained from UCSCs Table Browser utility under "Regulation" group and saving the output as a bed file. You can also obtain it from UCSCs database under the "cpgIslandExt.txt.gz" file, but be sure to convert the .txt file into a proper .bed file before usage.',
                    required=False)

parser.add_argument('-var', '--variation', dest="var_files", default=False,
                    help="Annotation file containing variation regions found in the genome (can be SNPs, strucutural variations, mutations or custom annotations). Can be obtained from UCSCs database or from Ensembl's GVF ftp directory.  Can take multiple files as input by using wildcard characters (*) with single quotation marks (ex. 'Homo_sapiens_*.gvf.gz'). Default: False",
                    required=False)

parser.add_argument('-ese', '--exon_splicing_enhancer', dest="ese_list",
                    default=False,
                    help="Txt list containing k-mers for exon splicing enhancer sequences. Default: False",
                    metavar="ese.txt", required=False)

parser.add_argument('-ess', '--exon_splicing_silencer', dest="ess_list",
                    default=False,
                    help="Txt list containing k-mers for exon splicing silencer sequences. Default: False",
                    metavar="ess.txt", required=False)

parser.add_argument('-ise', '--intron_splicing_enhancer', dest="ise_list",
                    default=False,
                    help="Txt list containing k-mers for intron splicing enhancer sequences. Default: False",
                    metavar="ise.txt", required=False)

parser.add_argument('-iss', '--intron_splicing_silencer', dest="iss_list",
                    default=False,
                    help="Txt list containing k-mers for intron splicing silencer sequences. Default: False",
                    metavar="iss.txt", required=False)

parser.add_argument('-custom_exon', '--custom-sequence-exon',
                    dest="custom_seq_ex", default=False,
                    help="Txt list containing custom sequences to be counted in the exons. Default: False",
                    metavar="custom_exon.txt", required=False)

parser.add_argument('-custom_intron', '--custom-sequence-intron',
                    dest="custom_seq_int", default=False,
                    help="Txt list containing custom sequences to be counted in the introns. Default: False",
                    metavar="custom_intron.txt", required=False)

parser.add_argument('-k', '--kmer', nargs='+', dest="kmer_list",
                    help="List of INT to create k-mers for countin. Default: False",
                    type=int, default=False, required=False)

parser.add_argument('-nuccont', '--nucleotide_content', dest="nuc_info",
                    default=2, metavar="nuc", required=False,
                    help="Defines the ammount of information included from the nucleotide sequence, 3 options available: 'Simple','Intermediate','Full'. Options:1 = Simple:[Length and pGC], 2 = Intermediate:[Length, pGC, pG, pC, pA, pT], 3 = Full:[All data from BedTools nucleotide sequence].' Default: 2 (Intermediate); p = percentage")

parser.add_argument("--maxEntScan", dest="max_ent_scan",
                    action="store_true", default=True,
                    help="Use this option if you want to use maxEntScan algorithm for scoring 3' and 5' splice sites. Although the splice sites are generally conserved, the original algorithm is trained for human datasets, which can lead to shaky results for other species. Default: True")

parser.add_argument("--keepBED", dest="keep_bed",
                    action="store_true", default=False,
                    help="Save the bed files generated for each class along with their up/downstream files. Default: False")

parser.add_argument("--keepTEMP", dest="keep_temp",
                    action="store_true", default=False,
                    help="Keep the temporary files generated for each transcript during up/downstream search. Default: False")

parser.add_argument("--ncores", dest="ncores", default=(mp.cpu_count() - 1),
                    help="Number of CPU cores used to process downtream/upstream exon search. Default:(ALL_CORES)-1",
                    type=int, metavar='INT')

parser.add_argument("--useThreads", dest="use_threads",
                    action="store_true", default=False,
                    help="Spawns multiple threads instead of multiple processes for parallelization. Is slower and uses mor RAM, but requires less CPU power. Default: False")

parser.add_argument("--debug", dest="debug", metavar='INT',
                    type=int, default=False,
                    help="Only seaches for the first N transcripts. Useful for debugin and checking if the code is working properly. Default: False")

args = parser.parse_args()


## Load the reference GTF file and convert it to a DataFrame

if not args.debug:
    pass
else:
    print()
    print("Running in debug mode. Only the first " + str(
        args.debug) + " entries will be used.")

if not args.use_threads:
    pass
else:
    print()
    print("Running with multiple threads. Watch out for that RAM!")

# if args.genome_browser == False:
#    pass
# else:
#    print
#    print "Making bed output compatible with UCSC's Genome Browser (WARNING: Remember to use UCSC's FASTA file instead of #ENSEMBL's due to differences in chromosome annotation)."


def get_transcript_id(df):
    a = df.copy()
    a['transcript_id'] = a['name'].apply(
        lambda x: x.split('transcript_id_')[1])
    a.sort_values(by=['transcript_id'], inplace=True)
    a.reset_index(drop=True, inplace=True)
    return a


def sorted_search(df, pattern):
    a = int(df['transcript_id'].searchsorted(pattern, 'left'))
    b = int(df['transcript_id'].searchsorted(pattern, 'right'))
    cut = df.iloc[a:b]
    cut.drop('transcript_id', 1, inplace=True)
    return cut

def get_upstream_and_downstream_exons(transcript_id):
    ##Select the exons of the transcript based on the list created above, then use this information
    ##to get some data of the transcript for further analysis and convert it into a BedTool object

    transcript_all_exons = sorted_search(all_exons_df, transcript_id)

    exon_count = transcript_all_exons.shape[0]  # number of exons
    strand_counts = transcript_all_exons['strand'].value_counts().shape[
        0]  # how many strands
    strand_id = \
        pd.DataFrame(transcript_all_exons['strand']).reset_index()['strand'][
            0]  # which strand

    all_exons_bedtool = BedTool.from_dataframe(transcript_all_exons).sort()

    if (strand_counts == 1) and (strand_id == "+"):
        if exon_count == 2:
            last_exons_bedtool = BedTool.from_dataframe(
                sorted_search(last_exons_df, transcript_id)).sort()
            first_exons_bedtool = BedTool.from_dataframe(
                sorted_search(first_exons_df, transcript_id)).sort()

            # TODO: refactor common closest call
            last_exons_bedtool.closest(all_exons_bedtool, s=True, id=True,
                                       D='ref', io=True, N=False,
                                       header=True,
                                       output=transcript_id + '_last_exons_up_temp.table')

            # TODO: refactor common closest call
            first_exons_bedtool.closest(all_exons_bedtool, s=True, iu=True,
                                        D='ref', io=True, N=False,
                                        header=True,
                                        output=transcript_id + '_first_exons_dn_temp.table')

        if exon_count >= 3:
            last_exons_bedtool = BedTool.from_dataframe(
                sorted_search(last_exons_df, transcript_id)).sort()
            first_exons_bedtool = BedTool.from_dataframe(
                sorted_search(first_exons_df, transcript_id)).sort()
            middle_exons_bedtool = BedTool.from_dataframe(
                sorted_search(middle_exons_df, transcript_id)).sort()

            # TODO: refactor common closest call
            last_exons_bedtool.closest(all_exons_bedtool,
                                       s=True,
                                       id=True,
                                       D='ref',
                                       io=True,
                                       N=False,
                                       header=True,
                                       output=transcript_id + '_last_exons_up_temp.table')

            # TODO: refactor common closest call
            first_exons_bedtool.closest(all_exons_bedtool,
                                        s=True,
                                        iu=True,
                                        D='ref',
                                        io=True,
                                        N=False,
                                        header=True,
                                        output=transcript_id + '_first_exons_dn_temp.table')

            # TODO: refactor common closest call
            middle_exons_bedtool.closest(all_exons_bedtool,
                                         s=True,
                                         id=True,
                                         D='ref',
                                         io=True,
                                         N=False,
                                         header=True,
                                         output=transcript_id + '_middle_exons_up_temp.table')

            # TODO: refactor common closest call
            middle_exons_bedtool.closest(all_exons_bedtool,
                                         s=True,
                                         iu=True,
                                         D='ref',
                                         io=True,
                                         N=False,
                                         header=True,
                                         output=transcript_id + '_middle_exons_dn_temp.table')

    if (strand_counts == 1) and (strand_id == "-"):
        if exon_count == 2:
            last_exons_bedtool = BedTool.from_dataframe(
                sorted_search(last_exons_df, transcript_id)).sort()
            first_exons_bedtool = BedTool.from_dataframe(
                sorted_search(first_exons_df, transcript_id)).sort()

            last_exons_bedtool.closest(all_exons_bedtool, s=True, iu=True,
                                       D='ref', io=True, N=False,
                                       header=True,
                                       output=transcript_id + '_last_exons_up_temp.table')

            first_exons_bedtool.closest(all_exons_bedtool, s=True, id=True,
                                        D='ref', io=True, N=False,
                                        header=True,
                                        output=transcript_id + '_first_exons_dn_temp.table')

        if exon_count >= 3:
            last_exons_bedtool = BedTool.from_dataframe(
                sorted_search(last_exons_df, transcript_id)).sort()
            first_exons_bedtool = BedTool.from_dataframe(
                sorted_search(first_exons_df, transcript_id)).sort()
            middle_exons_bedtool = BedTool.from_dataframe(
                sorted_search(middle_exons_df, transcript_id)).sort()

            last_exons_bedtool.closest(all_exons_bedtool,
                                       s=True,
                                       iu=True,
                                       D='ref',
                                       io=True,
                                       N=False,
                                       header=True,
                                       output=transcript_id + '_last_exons_up_temp.table')

            first_exons_bedtool.closest(all_exons_bedtool,
                                        s=True,
                                        id=True,
                                        D='ref',
                                        io=True,
                                        N=False,
                                        header=True,
                                        output=transcript_id + '_first_exons_dn_temp.table')

            middle_exons_bedtool.closest(all_exons_bedtool,
                                         s=True,
                                         iu=True,
                                         D='ref',
                                         io=True,
                                         N=False,
                                         header=True,
                                         output=transcript_id + '_middle_exons_up_temp.table')

            middle_exons_bedtool.closest(all_exons_bedtool,
                                         s=True,
                                         id=True,
                                         D='ref',
                                         io=True,
                                         N=False,
                                         header=True,
                                         output=transcript_id + '_middle_exons_dn_temp.table')

    if strand_counts == 2:
        print(
            'Both strands are present in transcript ' + transcript_id + ', check if all exons are in the same strand.')

def concatFiles(input_files, output_file):
    files = glob.glob(input_files)
    with open(output_file, "wb") as fo:
        for f in sorted(files):
            with open(f, "rb") as fi:
                shutil.copyfileobj(fi, fo)
                
def get_upstream_features(table):
    cat = pd.read_table(table, header=None)
    cat['tag'] = cat[
        3]  # +'_'+cat[9] ##The second part concatenates the upstream exon ID, but for simplicity we will not use that portion. If you want this information just uncoment the line and make the necessary changes to the script.

    cat_p = cat[cat[5] == '+'].copy()

    df_introns_p = cat_p[[0, 8, 1, 'tag', 4, 5]]
    df_introns_p[8] = cat_p[8]
    df_introns_p[1] = cat_p[1]
    df_introns_p = df_introns_p.rename(
        columns={0: 'chrom', 8: 'start', 1: 'end',
                 'tag': 'name', 4: 'score', 5: 'strand'})  # .drop_duplicates()

    ss_5p_p = cat_p[[0, 8, 1, 'tag', 4, 5]]
    ss_5p_p[8] = cat_p[8] - 3
    ss_5p_p[1] = cat_p[8] + 6
    ss_5p_p = ss_5p_p.rename(columns={0: 'chrom', 8: 'start', 1: 'end',
                                      'tag': 'name', 4: 'score', 5: 'strand'})

    ss_3p_p = cat_p[[0, 8, 1, 'tag', 4, 5]]
    ss_3p_p[8] = cat_p[1] - 20
    ss_3p_p[1] = cat_p[1] + 3
    ss_3p_p = ss_3p_p.rename(columns={0: 'chrom', 8: 'start', 1: 'end',
                                      'tag': 'name', 4: 'score', 5: 'strand'})

    bp_seq_p = cat_p[[0, 8, 1, 'tag', 4, 5]]
    bp_seq_p['length'] = bp_seq_p[1] - bp_seq_p[8]
    bp_seq_p_200 = bp_seq_p[bp_seq_p['length'] >= 200]
    bp_seq_p_small = bp_seq_p[bp_seq_p['length'] < 200]

    bp_seq_p_200[8] = cat_p[1] - 101
    bp_seq_p_200[1] = cat_p[1]

    bp_seq_p_small[8] = cat_p[1] - bp_seq_p['length']
    bp_seq_p_small[1] = cat_p[1]

    bp_seq_p = pd.concat([bp_seq_p_200, bp_seq_p_small]
                         ).reindex(columns=bp_seq_p_200.columns
                                   ).rename(
        columns={0: 'chrom', 8: 'start', 1: 'end',
                 'tag': 'name', 4: 'score', 5: 'strand'})

    cat_m = cat[cat[5] == '-'].copy()

    df_introns_m = cat_m[[0, 2, 7, 'tag', 4, 5]]
    df_introns_m[2] = cat_m[2]
    df_introns_m[7] = cat_m[7]
    df_introns_m = df_introns_m.rename(
        columns={0: 'chrom', 2: 'start', 7: 'end',
                 'tag': 'name', 4: 'score', 5: 'strand'})  # .drop_duplicates()

    ss_5p_m = cat_m[[0, 2, 7, 'tag', 4, 5]]
    ss_5p_m[2] = cat_m[7] - 6
    ss_5p_m[7] = cat_m[7] + 3
    ss_5p_m = ss_5p_m.rename(columns={0: 'chrom', 'tag': 'name', 2: 'start',
                                      7: 'end', 4: 'score', 5: 'strand'})

    ss_3p_m = cat_m[[0, 2, 7, 'tag', 4, 5]]
    ss_3p_m[2] = cat_m[2] - 3
    ss_3p_m[7] = cat_m[2] + 20
    ss_3p_m = ss_3p_m.rename(columns={0: 'chrom', 'tag': 'name', 2: 'start',
                                      7: 'end', 4: 'score', 5: 'strand'})

    bp_seq_m = cat_m[[0, 2, 7, 'tag', 4, 5]].copy()
    bp_seq_m['length'] = cat_m[7] - cat_m[2]
    bp_seq_m_200 = bp_seq_m[bp_seq_m['length'] >= 200]
    bp_seq_m_small = bp_seq_m[bp_seq_m['length'] < 200]

    bp_seq_m_200[2] = cat_m[2]
    bp_seq_m_200[7] = cat_m[2] + 100

    bp_seq_m_small[2] = cat_m[2]
    bp_seq_m_small[7] = cat_m[2] + bp_seq_m['length']

    bp_seq_m = pd.concat([bp_seq_m_200, bp_seq_m_small]
                         ).reindex(columns=bp_seq_m_200.columns
                                   ).rename(
        columns={0: 'chrom', 2: 'start', 7: 'end',
                 'tag': 'name', 4: 'score', 5: 'strand'})

    df_exons = cat[[6, 7, 8, 'tag', 10, 11]].drop_duplicates().rename(
        columns={6: 'chrom', 7: 'start', 8: 'end', 'tag': 'name', 10: 'score',
                 11: 'strand'})

    df_introns = pd.concat([df_introns_p, df_introns_m]
                           ).drop_duplicates().reindex(
        columns=df_introns_p.columns).dropna()

    df_5p_ss = pd.concat([ss_5p_p, ss_5p_m]
                         ).drop_duplicates().reindex(
        columns=ss_5p_p.columns).dropna()

    df_3p_ss = pd.concat([ss_3p_p, ss_3p_m]
                         ).drop_duplicates().reindex(
        columns=ss_3p_p.columns).dropna()

    df_bp_seq = pd.concat([bp_seq_p, bp_seq_m]
                          ).drop_duplicates().reindex(
        columns=bp_seq_p.columns).drop('length', 1).dropna()

    return df_exons, df_introns, df_5p_ss, df_3p_ss, df_bp_seq

def get_downstream_features(table):
    cat = pd.read_table(table, header=None)
    cat['tag'] = cat[
        3]  # +'_'+cat[9] ##The second part concatenates the upstream exon ID, but for simplicity we will not use that portion. If you want this information just uncoment the line and make the necessary changes to the script.

    df_exons = cat[[6, 7, 8, 'tag', 10, 11]].drop_duplicates().rename(
        columns={6: 'chrom', 7: 'start', 8: 'end', 'tag': 'name', 10: 'score',
                 11: 'strand'})

    cat_p = cat[cat[5] == '+'].copy()

    df_introns_p = cat_p[[0, 2, 7, 'tag', 4, 5]]
    df_introns_p[2] = cat_p[2]
    df_introns_p[7] = cat_p[7]
    df_introns_p = df_introns_p.rename(
        columns={0: 'chrom', 2: 'start', 7: 'end',
                 'tag': 'name', 4: 'score', 5: 'strand'})

    ss_5p_p = cat_p[[0, 2, 7, 'tag', 4, 5]]
    ss_5p_p[2] = cat_p[2] - 3
    ss_5p_p[7] = cat_p[2] + 6
    ss_5p_p = ss_5p_p.rename(columns={0: 'chrom', 2: 'start', 7: 'end',
                                      'tag': 'name', 4: 'score', 5: 'strand'})

    ss_3p_p = cat_p[[0, 2, 7, 'tag', 4, 5]]
    ss_3p_p[2] = cat_p[7] - 20
    ss_3p_p[7] = cat_p[7] + 3
    ss_3p_p = ss_3p_p.rename(columns={0: 'chrom', 2: 'start', 7: 'end',
                                      'tag': 'name', 4: 'score', 5: 'strand'})

    bp_seq_p = cat_p[[0, 2, 7, 'tag', 4, 5]]
    bp_seq_p['length'] = bp_seq_p[7] - bp_seq_p[2]
    bp_seq_p_200 = bp_seq_p[bp_seq_p['length'] >= 200]
    bp_seq_p_small = bp_seq_p[bp_seq_p['length'] < 200]

    bp_seq_p_200[2] = cat_p[7] - 100
    bp_seq_p_200[7] = cat_p[7]

    bp_seq_p_small[2] = cat_p[7] - bp_seq_p['length']
    bp_seq_p_small[7] = cat_p[7]

    bp_seq_p = pd.concat([bp_seq_p_200, bp_seq_p_small]
                         ).reindex(columns=bp_seq_p_200.columns
                                   ).rename(
        columns={0: 'chrom', 2: 'start', 7: 'end',
                 'tag': 'name', 4: 'score', 5: 'strand'})

    cat_m = cat[cat[5] == '-'].copy()
    df_introns_m = cat_m[[0, 8, 1, 'tag', 4, 5]]
    df_introns_m[8] = cat_m[8]
    df_introns_m[1] = cat_m[1]
    df_introns_m = df_introns_m.rename(
        columns={0: 'chrom', 8: 'start', 1: 'end',
                 'tag': 'name', 4: 'score', 5: 'strand'})

    ss_5p_m = cat_m[[0, 8, 1, 'tag', 4, 5]]
    ss_5p_m[8] = cat_m[1] - 6
    ss_5p_m[1] = cat_m[1] + 3
    ss_5p_m = ss_5p_m.rename(columns={0: 'chrom', 'tag': 'name', 8: 'start',
                                      1: 'end', 4: 'score', 5: 'strand'})

    ss_3p_m = cat_m[[0, 8, 1, 'tag', 4, 5]]
    ss_3p_m[8] = cat_m[8] - 3
    ss_3p_m[1] = cat_m[8] + 20
    ss_3p_m = ss_3p_m.rename(columns={0: 'chrom', 'tag': 'name', 8: 'start',
                                      1: 'end', 4: 'score', 5: 'strand'})

    bp_seq_m = cat_m[[0, 8, 1, 'tag', 4, 5]].copy()
    bp_seq_m['length'] = cat_m[1] - cat_m[8]
    bp_seq_m_200 = bp_seq_m[bp_seq_m['length'] >= 200]
    bp_seq_m_small = bp_seq_m[bp_seq_m['length'] < 200]

    bp_seq_m_200[8] = cat_m[8]
    bp_seq_m_200[1] = cat_m[8] + 100

    bp_seq_m_small[8] = cat_m[8]
    bp_seq_m_small[1] = cat_m[8] + bp_seq_m['length']

    bp_seq_m = pd.concat([bp_seq_m_200, bp_seq_m_small]
                         ).reindex(columns=bp_seq_m_200.columns
                                   ).rename(
        columns={0: 'chrom', 8: 'start', 1: 'end',
                 'tag': 'name', 4: 'score', 5: 'strand'})

    df_introns = pd.concat([df_introns_p, df_introns_m]
                           ).drop_duplicates().reindex(
        columns=df_introns_p.columns).dropna()

    df_5p_ss = pd.concat([ss_5p_p, ss_5p_m]
                         ).drop_duplicates().reindex(
        columns=ss_5p_p.columns).dropna()

    df_3p_ss = pd.concat([ss_3p_p, ss_3p_m]
                         ).drop_duplicates().reindex(
        columns=ss_3p_p.columns).dropna()

    df_bp_seq = pd.concat([bp_seq_p, bp_seq_m]
                          ).drop_duplicates().reindex(
        columns=bp_seq_p.columns).drop('length', 1).dropna()

    return df_exons, df_introns, df_5p_ss, df_3p_ss, df_bp_seq

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

def splice_site_fasta(df, name):
    bt = BedTool.from_dataframe(df)  # .remove_invalid().saveas()
    bt.sequence(fi=genome_fasta, fo=name + ".fasta", name=True, s=True)
    return name + ".fasta"


def run_maxent3p(fasta, cmd="perl maxEntScan/score3.pl"):
    composed_command = " ".join([cmd, fasta])
    p = Popen(composed_command, stdout=PIPE, shell=True)
    stdout, stderr = p.communicate()
    return pd.read_table(StringIO(stdout), sep='\t', header=None
                         ).rename(columns={1: 'name', 2: 'ss_scr'}).drop(0, 1)

def run_maxent5p(fasta, cmd="perl maxEntScan/score5.pl"):
    composed_command = " ".join([cmd, fasta])
    p = Popen(composed_command, stdout=PIPE, shell=True)
    stdout, stderr = p.communicate()
    return pd.read_table(StringIO(stdout), sep='\t', header=None
                         ).rename(columns={1: 'name', 2: 'ss_scr'}).drop(0, 1)

def get_maxent_score(df, name):
    if name.find("_3ss") != -1:
        result = run_maxent3p(splice_site_fasta(df, name))
        result['name'] = result['name'].apply(lambda x: x.split('::')[0])
        return result
    if name.find("_5ss") != -1:
        result = run_maxent5p(splice_site_fasta(df, name))
        result['name'] = result['name'].apply(lambda x: x.split('::')[0])
        return result
    else:
        pass

def get_conservation_scores(con_file, df, cmd="bigWigAverageOverBed"):
    print()
    print("Getting conservation from: " + str(con_file))
    print()
    BedTool.from_dataframe(df).saveas('conservation_temp.bed')
    composed_command = " ".join([cmd, con_file, 'conservation_temp.bed',
                                 'conservation_result_temp.tab'])
    p = Popen(composed_command, shell=True)
    stdout, stderr = p.communicate()
    source = str(con_file.split('.')[1])
    result = pd.concat(pd.read_table('conservation_result_temp.tab',
                                     names=('name', 'length', 'covered',
                                            'sum', 'mean_' + source,
                                            'mean0_' + source), iterator=True,
                                     chunksize=10000
                                     ), ignore_index=True)
    os.remove('conservation_temp.bed')
    os.remove('conservation_result_temp.tab')
    return result[['name', 'mean_' + source]]

def get_cpg_islands(bedtool):
    cpg = BedTool(args.cpg_file)
    cpg_counts = pd.concat(
        bedtool.intersect(cpg, c=True).to_dataframe(iterator=True,
                                                    chunksize=10000
                                                    ),
        ignore_index=True).rename(columns={'thickStart': 'CpG_count'})
    return cpg_counts

def get_var_counts(bedtool, var_file):
    print()
    print("Getting variation from: " + str(var_file))
    print()
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

def get_kmer_counts(df):
    print
    print "Extracting K-mer counts"
    bases=['A','T','G','C']
    kmers = pd.DataFrame()
    for i in range(len(args.kmer_list)):
        k = pd.DataFrame([''.join(p) for p in itertools.product(bases, repeat=args.kmer_list[i])])
        kmers = pd.concat([kmers, k])
    kmers.reset_index(inplace=True)
    kmers.drop('index', 1, inplace=True)
    for i in range(len(kmers[0])):
        df['kmer_'+str(kmers[0][i])] = df.apply(lambda x: str(x['seq']).upper().count(str(kmers[0][i]).upper()),1)
    return df

def get_data(df, name, matrix):
    print()
    print("Starting " + name)

    bedtool = BedTool.from_dataframe(df).saveas()
    a = nuc_cont(bedtool)
    
    if args.kmer_list:
        a = get_kmer_counts(a)
    elif not args.kmer_list:
        pass

    if args.cpg_file:
        cpg_counts = get_cpg_islands(bedtool)
        a = a.merge(cpg_counts[['name', 'CpG_count']],
                    on='name').drop_duplicates()
    elif not args.cpg_file:
        pass

    if args.var_files:
        var_files_list = glob.glob(args.var_files)
        for i in range(len(var_files_list)):
            var = get_var_counts(bedtool, var_files_list[i])
            a = a.merge(var, on='name')
    elif not args.var_files:
        pass

    if args.con_files:
        con_files_list = glob.glob(args.con_files)
        for i in range(len(con_files_list)):
            con = get_conservation_scores(con_files_list[i], df)
            a = a.merge(con, on='name')
    elif not args.con_files:
        pass

    if name.find("_exon") != -1 and name.find('single') == -1:
        if not args.ese_list:
            pass
        elif args.ese_list:
            a['ese_count'] = a['seq'].apply(
                lambda x: sum(x.upper().count(y.upper()) for y in ese_list))

        if not args.ess_list:
            pass
        elif args.ess_list:
            a['ess_count'] = a['seq'].apply(
                lambda x: sum(x.upper().count(y.upper()) for y in ess_list))

        if not args.custom_seq_ex:
            pass
        elif args.custom_seq_ex:
            a['custom_seq_exon'] = a['seq'].apply(lambda x: sum(
                x.upper().count(y.upper()) for y in custom_seq_ex))

    if name.find("_intron") != -1:
        if not args.ise_list:
            pass
        elif args.ise_list:
            a['ise_count'] = a['seq'].apply(
                lambda x: sum(x.upper().count(y.upper()) for y in ise_list))

        if not args.iss_list:
            pass
        elif args.iss_list:
            a['iss_count'] = a['seq'].apply(
                lambda x: sum(x.upper().count(y.upper()) for y in iss_list))

        if not args.custom_seq_int:
            pass
        elif args.custom_seq_int:
            a['custom_seq_intron'] = a['seq'].apply(lambda x: sum(
                x.upper().count(y.upper()) for y in custom_seq_int))

    if args.max_ent_scan:
        if name.find("_3ss") != -1 or name.find("_5ss") != -1:
            ss = get_maxent_score(df, name)
            a = a.merge(ss, on='name').drop_duplicates()
        else:
            pass
    else:
        pass

    if str(args.nuc_info) == 3:
        z = a.drop(['seqname', 'start', 'end', 'score', 'strand', 'seq'], 1)
    else:
        z = a.drop(['seq'], 1)

    z = z.set_index('name').add_suffix('_' + name).reset_index()
    z = filter_columns(z)
    print()
    print(name + " finished")
    return matrix.merge(z, on='name').drop_duplicates()    

##Load the genome file that matches the version of the GTF you are using. Pysam will be used to build an index of
##the FASTA file.

genome_fasta = args.genome_file
pysam.faidx(genome_fasta)

## Load the GTF file and convert it to a DataFrame

print()
print("Loading GTF file and creating the DataFrame")

gtf = BedTool(args.gtf_file)
df = pd.concat(gtf.to_dataframe(iterator=True, chunksize=10000
                                ), ignore_index=True).dropna().sort_values(
    by=['seqname', 'start'])

##Subtract 1 from the start coordinate of every sequence to facilitate the conversion to bed format later.
##This is because bed is 0-based while gtf is 1-based (IMPORTANT!!!!)
##Then extract the lines which contain transcripts and exons, add a "tag" annotation
##Warning: this is designed to work on GTF files downloaded from ensembl. If your GTF file already contains 
##the 'chr' annotation, just comment out the line "df['seqname'] = 'chr'+df['seqname'].astype(str)"
##Warning2: If you would like to export the bed files generated by the code, add the chr annotation and
##remove the nonchromossomal annotation from the file. To do this, just uncomment the lines 335-337 in the code

#    df = df[~df['seqname'].str.contains('GL|KI').fillna(False)]
#    df['seqname'] = 'chr'+df['seqname'].astype(str)
#    df['seqname'] = df['seqname'].str.replace('chrMT','chrM')

df['start'] = (df['start'] - 1).astype(int)
df['end'] = df['end'].astype(int)

print()
print("Filtering exons annotations")

df_exons = df[df['feature'] == 'exon'].reset_index().drop('index', 1)
df_exons['exon_id'] = 'exon_id_E0' + (df_exons.index + 1).astype(str)
df_exons['transcript_id'] = "transcript_id_" + df_exons['attributes'].apply(
    lambda x: x.split('transcript_id "')[1])
df_exons['transcript_id'] = df_exons['transcript_id'].apply(
    lambda x: x.split('"')[0])

df_exons['tag'] = df_exons['exon_id'] + "_" + df_exons['transcript_id']
df_exons['score'] = 0

print()
print(
    "Classifing exons in: single exons, first exons, last exons and middle exons.")

## Select only the exons which match the start (first exons) and end (last exons) of the transcript and convert them to
## BED format.

## Single exons
exons_grouped = df_exons.groupby('transcript_id')
group_size_filter = (exons_grouped.size() == 1)
single_exons = df_exons[
    df_exons['transcript_id'].isin(group_size_filter[group_size_filter].index)]

single_exons_bed = single_exons[
    ['seqname', 'start', 'end', 'tag', 'score', 'strand']
].drop_duplicates().sort_values(by=['seqname', 'start']
                                ).rename(columns={'tag': 'name'})

##Create a frame for non-single exons:

df_non_single_exons = df_exons[~df_exons['tag'].isin(single_exons['tag'])]

## For first and last exons we need to separate exons by strand

df_ns_exons_plus = df_non_single_exons[df_non_single_exons['strand'] == '+']
df_ns_exons_minus = df_non_single_exons[df_non_single_exons['strand'] == '-']

# First Exons

first_exons_plus = df_ns_exons_plus.groupby(
    ['transcript_id']).first().reset_index()
first_exons_minus = df_ns_exons_minus.groupby(
    ['transcript_id']).last().reset_index()

first_exons = pd.concat([first_exons_plus, first_exons_minus]).reindex(
    columns=df_exons.columns)

first_exons_bed = first_exons[
    ['seqname', 'start', 'end', 'tag', 'score', 'strand']
].drop_duplicates().sort_values(by=['seqname', 'start']
                                ).rename(columns={'tag': 'name'})

# Last Exons

last_exons_plus = df_ns_exons_plus.groupby(
    ['transcript_id']).last().reset_index()
last_exons_minus = df_ns_exons_minus.groupby(
    ['transcript_id']).first().reset_index()

last_exons = pd.concat([last_exons_plus, last_exons_minus]).reindex(
    columns=df_exons.columns)

last_exons_bed = last_exons[
    ['seqname', 'start', 'end', 'tag', 'score', 'strand']
].drop_duplicates().sort_values(by=['seqname', 'start']
                                ).rename(columns={'tag': 'name'})

## Get middle exons by removing the first and last exons from the total found, then create a bed file.

middle_exons = df_non_single_exons[
    (~(df_non_single_exons['tag'].isin(first_exons['tag'])) &
     ~(df_non_single_exons['tag'].isin(last_exons['tag'])) &
     ~(df_non_single_exons['tag'].isin(single_exons['tag'])))]

middle_exons_bed = middle_exons[
    ['seqname', 'start', 'end', 'tag', 'score', 'strand']
].drop_duplicates().sort_values(by=['seqname', 'start']
                                ).rename(columns={'tag': 'name'})


## Apply the same process for creating a BED file for all exons

all_exons_bed = df_exons[['seqname', 'start', 'end', 'tag', 'score', 'strand']
].drop_duplicates().sort_values(by=['seqname', 'start']
                                ).rename(columns={'tag': 'name'})


## Delete the unsorted bedfiles and cleanup temporary pybedtool objects

pybedtools.helpers.cleanup()


## Now, we find which are the upstream and downstream exons (in the mRNA) for the exons before. For middle
## exons we select both upstream and downstream exons. using the "closest" tool from BedTool suite first
## exons only the downstream and for last exons only the upstream. A simples schematic of selection follows:
##
## Middle exons (+ strand): E1 (up_exon) ----- E2 (middle_exon) ----- E3 (dn_exon)
## Middle exons (- strand): E1 (dn_exon) ----- E2 (middle_exon) ----- E3 (up_exon)
##
## First exons (+ strand): E1 (first_exon) ----- E2 (dn_exon)
## First exons (- strand): E2 (dn_exon) ----- E3 (first_exon)
##
## Last exons (+ strand): E2 (up_exon) ----- E3 (last_exon)
## Last exons (- strand): E1 (last_exon) ----- E2 (up_exon)
##
## Addiotionally, by using the relative coordinates of the exons, we can infer the intronic coordinates. 
## Using the start/end coordinates of the exons and their upstream/downstream partners, adding +1 and -1 
## for correction. For 5' splice sites scoring by MaxEntScan, we need to get 9nt FASTA files that partially 
## inside the exon (3nt) and part in the intron (6 nt) (ex CAGgtaagt). For analysis of the 3' splice site, 
## we also need to get a hybrid FASTA file (20nt in intron and 3nt in exon) (ex. ttccaaacgaacttttgtagGGA). 
## For branch-point prediction using SVM-BP finder, we need the 200 nt on the 3' end of the intron. 
## To obtain the necessary coordinates, we will manipulate the coordinates obtained from the upstream 
## and downstream exons.

## To get, for each transcript, it's upstream and downstream exon, we nedd to create a function. 
## Since this would take a long time we need to parallellize in multiple cores. 
## First we find which are the non-single exon transcripts:

all_exons_df = get_transcript_id(all_exons_bed)
first_exons_df = get_transcript_id(first_exons_bed)
last_exons_df = get_transcript_id(last_exons_bed)
middle_exons_df = get_transcript_id(middle_exons_bed)
single_exons_df = get_transcript_id(single_exons_bed)

all_transcripts = all_exons_df['transcript_id']
single_transcripts = single_exons_bed['name'].apply(
    lambda x: x.split('transcript_id_')[1])

if not args.debug:
    non_single_transcripts = all_transcripts[
        all_transcripts.isin(single_transcripts) == False
        ].drop_duplicates().reset_index(drop=True)
else:
    non_single_transcripts = all_transcripts[
        all_transcripts.isin(single_transcripts) == False
        ].drop_duplicates().reset_index(drop=True).head(args.debug)


# Now we use a function that, for each transcript, find all the exons of that and run bedtool closest
# To get it's upstream and downstream neighbours

print()
print(
    "Searching for upstream and downstream exons in each transcript. (WARNING: This may take a while, go grab a cup of cofee... or a movie... or go home and get some sleep, it'll take a few hours to run)")

if args.use_threads:
    if __name__ == '__main__':
        p = ThreadPool(args.ncores)
        p.map(get_upstream_and_downstream_exons, non_single_transcripts)

elif not args.use_threads:
    if __name__ == '__main__':
        p = Pool(args.ncores)
        p.map(get_upstream_and_downstream_exons, non_single_transcripts)

concatFiles('*_last_exons_up_temp.table', 'upstream_last_exons_temp.table')
concatFiles('*_first_exons_dn_temp.table', 'downstream_first_exons_temp.table')
concatFiles('*_middle_exons_up_temp.table', 'upstream_middle_exons_temp.table')
concatFiles('*_middle_exons_dn_temp.table',
            'downstream_middle_exons_temp.table')


## Find the regions

## Find up/down exons for Middle exons

print()
print('Finding upstream features for middle exons')

up_middle_exon, up_middle_intron, up_middle_5ss, up_middle_3ss, up_middle_bp = get_upstream_features(
    'upstream_middle_exons_temp.table')

print()
print('Finding downstream features for middle exons')

dn_middle_exon, dn_middle_intron, dn_middle_5ss, dn_middle_3ss, dn_middle_bp = get_downstream_features(
    'downstream_middle_exons_temp.table')

## Find down exons for First exons

print()
print('Finding downstream features for first exons')

dn_first_exon, dn_first_intron, dn_first_5ss, dn_first_3ss, dn_first_bp = get_downstream_features(
    'downstream_first_exons_temp.table')

## Find up exons for Last exons

print()
print('Finding upstream features for last exons')

up_last_exon, up_last_intron, up_last_5ss, up_last_3ss, up_last_bp = get_upstream_features(
    'upstream_last_exons_temp.table')


## Save all of them to a .bed file

if args.keep_bed:
    print()
    print("Saving classified exons to bed files")
    save_bed(single_exons_bed, str(args.prefix) + '.single.exons.bed.gz')
    save_bed(first_exons_bed, str(args.prefix) + '.first.exons.bed.gz')
    save_bed(last_exons_bed, str(args.prefix) + '.last.exons.bed.gz')
    save_bed(middle_exons_bed, str(args.prefix) + '.middle.exons.bed.gz')
    save_bed(all_exons_bed, str(args.prefix) + '.all.exons.bed.gz')

    print()
    print("Saving features to bed files")

    # Up middle
    save_bed(up_middle_exon,
             str(args.prefix) + '.upstream_middle_exons.bed.gz')
    save_bed(up_middle_intron,
             str(args.prefix) + '.upstream_middle_introns.bed.gz')
    save_bed(up_middle_5ss, str(args.prefix) + '.upstream_middle_5ss.bed.gz')
    save_bed(up_middle_3ss, str(args.prefix) + '.upstream_middle_3ss.bed.gz')
    save_bed(up_middle_bp, str(args.prefix) + '.upstream_middle_bp.bed.gz')

    # Dn middle
    save_bed(dn_middle_exon,
             str(args.prefix) + '.downstream_middle_exons.bed.gz')
    save_bed(dn_middle_intron,
             str(args.prefix) + '.downstream_middle_introns.bed.gz')
    save_bed(dn_middle_5ss, str(args.prefix) + '.downstream_middle_5ss.bed.gz')
    save_bed(dn_middle_3ss, str(args.prefix) + '.downstream_middle_3ss.bed.gz')
    save_bed(dn_middle_bp, str(args.prefix) + '.downstream_middle_bp.bed.gz')

    # Dn first
    save_bed(dn_first_exon,
             str(args.prefix) + '.downstream_first_exons.bed.gz')
    save_bed(dn_first_intron,
             str(args.prefix) + '.downstream_first_introns.bed.gz')
    save_bed(dn_first_5ss, str(args.prefix) + '.downstream_first_5ss.bed.gz')
    save_bed(dn_first_3ss, str(args.prefix) + '.downstream_first_3ss.bed.gz')
    save_bed(dn_first_bp, str(args.prefix) + '.downstream_first_bp.bed.gz')

    # Up last
    save_bed(up_last_exon, str(args.prefix) + '.upstream_last_exons.bed.gz')
    save_bed(up_last_intron,
             str(args.prefix) + '.upstream_last_introns.bed.gz')
    save_bed(up_last_5ss, str(args.prefix) + '.upstream_last_5ss.bed.gz')
    save_bed(up_last_3ss, str(args.prefix) + '.upstream_last_3ss.bed.gz')
    save_bed(up_last_bp, str(args.prefix) + '.upstream_last_bp.bed.gz')

    archive = tarfile.open("bed_files.tar.gz", "w|gz")
    for filename in glob.glob('*.bed.gz'):
        archive.add(filename)
    archive.close()

    list(map(os.remove, glob.glob("*.bed.gz")))

elif not args.keep_bed:
    print()
    print("Saving bed file with annotated exons")
    save_bed(all_exons_bed, str(args.prefix) + '.all.exons.bed.gz')

if args.keep_temp:
    print()
    print("Compressing all temporary files and cleaning up")
    archive = tarfile.open("temp_tables.tar.gz", "w|gz")
    for filename in glob.glob('*.table'):
        archive.add(filename)
    archive.close()
    list(map(os.remove, glob.glob("*.table")))

if not args.keep_temp:
    print()
    print("Removing all temporary files")
    list(map(os.remove, glob.glob("*.table")))

##Now we will start gathering data of each region and assemble everything into a dataframe which will be passed to 
## the classifier in posterior steps. For this we will need to use several functions.

if not args.ese_list:
    pass
else:
    ese_list = list(pd.read_table(args.ese_list, header=None)[0])

if not args.ess_list:
    pass
else:
    ess_list = list(pd.read_table(args.ess_list, header=None)[0])

if not args.ise_list:
    pass
else:
    ise_list = list(pd.read_table(args.ise_list, header=None)[0])

if not args.iss_list:
    pass
else:
    iss_list = list(pd.read_table(args.iss_list, header=None)[0])

if not args.custom_seq_ex:
    pass
else:
    custom_seq_ex = list(pd.read_table(args.custom_seq_ex, header=None)[0])

if not args.custom_seq_int:
    pass
else:
    custom_seq_int = list(pd.read_table(args.custom_seq_int, header=None)[0])

##Now, we will create the spine for our matrix by selecting only the names of the exons used in each class.
##We will also need a list of the frames associated with that class and the names of the objects

if args.debug:
    first_exons_bed = first_exons_bed[
        first_exons_bed['name'].isin(dn_first_3ss['name'])]
    middle_exons_bed = middle_exons_bed[
        middle_exons_bed['name'].isin(dn_middle_3ss['name'])]
    last_exons_bed = last_exons_bed[
        last_exons_bed['name'].isin(up_last_3ss['name'])]
    single_exons_bed = single_exons_bed.head(args.debug)
elif not args.debug:
    pass

matrix_first = pd.DataFrame(first_exons_bed['name'])
frames_first = list((first_exons_bed,
                     dn_first_exon, dn_first_intron, dn_first_3ss,
                     dn_first_5ss))
names_first = list(('first_exons',
                    'dn_first_exon', 'dn_first_intron', 'dn_first_3ss',
                    'dn_first_5ss'))

matrix_middle = pd.DataFrame(middle_exons_bed['name'])
frames_middle = list((middle_exons_bed,
                      dn_middle_exon, dn_middle_intron, dn_middle_3ss,
                      dn_middle_5ss,
                      up_middle_exon, up_middle_intron, up_middle_3ss,
                      up_middle_5ss))
names_middle = list(('middle_exons',
                     'dn_middle_ex', 'dn_middle_int', 'dn_middle_3ss',
                     'dn_middle_5ss',
                     'up_middle_ex', 'up_middle_int', 'up_middle_3ss',
                     'up_middle_5ss'))

matrix_last = pd.DataFrame(last_exons_bed['name'])
frames_last = list((last_exons_bed,
                    up_last_exon, up_last_intron, up_last_3ss, up_last_5ss))
names_last = list(('last_exons',
                   'up_last_exon', 'up_last_intron', 'up_last_3ss',
                   'up_last_5ss'))

matrix_single = pd.DataFrame(single_exons_bed['name'])

print()
print("Starting matrix build.")
print()
print("Getting data for first exons.")

for i in range(len(frames_first)):
    matrix_first = get_data(frames_first[i], names_first[i], matrix_first)

matrix_first.set_index('name').drop_duplicates().to_csv(
    str(args.prefix) + '.first.exons.datamatrix.tsv', sep='\t')

print()
print("First exons DONE. Starting middle exons.")

for i in range(len(frames_middle)):
    matrix_middle = get_data(frames_middle[i], names_middle[i], matrix_middle)

matrix_middle.set_index('name').drop_duplicates().to_csv(
    str(args.prefix) + '.middle.exons.datamatrix.tsv', sep='\t')

print()
print("Middle exons DONE. Starting last exons.")

for i in range(len(frames_last)):
    matrix_last = get_data(frames_last[i], names_last[i], matrix_last)

matrix_last.set_index('name').drop_duplicates().to_csv(
    str(args.prefix) + '.last.exons.datamatrix.tsv', sep='\t')

print()
print("Last exons DONE. Starting single exons.")

matrix_single = get_data(single_exons_bed, 'single_exons', matrix_single)
matrix_single.set_index('name').drop_duplicates().to_csv(
    str(args.prefix) + '.single.exons.datamatrix.tsv', sep='\t')

print()
print('Single exons DONE.')

print()
print("Data matrices build complete")
