import argparse
import sys


def get_args():
    parse = argparse.ArgumentParser(description='Process args for rm_verify')
    parse.add_argument('-o', help='Output directory to write files to')
    parse.add_argument('-i', help='Input directory if files already \
                                  downloaded or for test data')
    parse.add_argument('-k', help='Path to kallisto index if already created')
    parse.add_argument(
        '-q', help='Path to coverted fastq files if already created')
    parse.add_argument(
        '-b', help='Path bowtie2 index if already created')
    parse = parse.parse_args()

    # output directory should be miniProject_Ethan_Holleman
    # make the directory is already does not exist

    return parse
