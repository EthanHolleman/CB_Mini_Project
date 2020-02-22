import subprocess
import os
import csv


def run_kallisto(k_index, query_file_paths, output_dir, 
                 kallisto_executable='kallisto', b=30, t=4):
    '''
    Given the path to a kallisto index and a list of query files
    takes the query filepaths uses subprocess to run analysis using 
    kallisto program. 
    '''
    # make new directory for each call? Have to test that out.
    kallisto_dirs = []
    for query_a, query_b in query_file_paths:
        output_dir = os.path.join(output_dir, os.path.basename(query_a) + '_kallisto')
        kallisto_dirs.append(output_dir)
        cmd = [kallisto_executable, 'quant', '-i', k_index, '-o', output_dir,
               '-b', b, '-t', t, query_a, query_b]
    
        cmd = [str(i) for i in cmd]  # convert everything to string
        subprocess.call(cmd)
    
    return kallisto_dirs
    
    # > kallisto quant -i index/index.idx   -o results/DRR002318 -b 30 -t 4  data/DRR002318_1.fastq.gz   data/DRR002318_2.fastq.gz
    # need to write some kind of log file that can pass onto sleuth

def make_kalisto_index(trans_file, output_dir, kallisto_executable='kallisto'):
    '''
    Given a transcriptome fasta file and an output dir runs kallisto via
    subprocess to create a new kallisto index from that fasta file. Returns the
    path to the newly created kallisto index as a string.

    TODO: Add way to access and edit the log file created when making the
    Kalisto index.

    '''
    output_file = os.path.join(output_dir, 'mp_kalisto_index.idx')
    cmd = [kallisto_executable, 'index', '-i', output_file, trans_file,
           '--make-unique']
    subprocess.call(cmd)
    return output_file