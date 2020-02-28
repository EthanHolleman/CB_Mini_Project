import argparse
import sys
import os

from data import if_not_dir_make


def get_args():
    parse = argparse.ArgumentParser(description='Process args for rm_verify')
    parse.add_argument('-o', default='./miniProject_Ethan_Holleman',
                       help='Output directory to write files to')
    parse.add_argument('-i', help='Input directory if files already \
                                  downloaded or for test data')
    parse.add_argument('-k', help='Path to kallisto index if already created')
    parse.add_argument(
        '-r', help='Path to kallisto ressults if already created')
    parse.add_argument(
        '-q', help='Path to coverted fastq files if already created')
    parse.add_argument(
        '-b', help='Path bowtie2 index if already created')
    parse.add_argument(
        '-s', help='Path to bowtie results if already ran')
    parse.add_argument(
        '-f', help='Path to big fasta file')
    parse.add_argument(
        '-a', help='Path to assembled genome if exists')
    parse.add_argument('-t', default=2, help='Number of threads')
    parse.add_argument('-l', help='path to write log file to')
    parse.add_argument('-local', default=0, help='If set to 1 runs BLAST search locally')
    parse.add_argument('-test', default=0, help='Run the program in with test data. You will still need to set -local 1 to run BLAST locally.')

    parse = parse.parse_args()
    cwd_name = if_not_dir_make(os.getcwd(), 'miniProject_Ethan_Holleman')
    if not parse.o:
        parse.o = cwd_name
    if not parse.l:
        parse.l = cwd_name
    if parse.test == 1 or parse.test == '1':
        parse.q = './test_data/SRA_to_FASTQ'
        
        

    # output directory should be miniProject_Ethan_Holleman

    return parse
