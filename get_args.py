import argparse
import sys

def get_args():
    parse = argparse.ArgumentParser(description='Process args for rm_verify')
    parse.add_argument('-o', help='Output directory to write files to')
    parse.add_argument('-i', help='Input directory if files already downloaded or for test data')

    parse = parse.parse_args()


    return parse
