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
    kallisto_dirs = []
    for query_a, query_b in query_file_paths:  # iterate by paired fastq
        results_dir = os.path.join(output_dir,  # result dir based on query a
                                   os.path.basename(query_a) + '_kallisto')
        kallisto_dirs.append(results_dir)  # store path
        cmd = [kallisto_executable, 'quant', '-i', k_index, '-o', results_dir,
               '-b', b, '-t', t, query_a, query_b]
    
        cmd = [str(i) for i in cmd]  # convert everything to string
        subprocess.call(cmd)  # run command
    return kallisto_dirs
    

def make_kalisto_index(trans_file, output_dir, kallisto_executable='kallisto',
                       kalliso_index_name='mp_kalisto_index.idx'):
    '''
    Given a transcriptome fasta file and an output dir runs kallisto via
    subprocess to create a new kallisto index from that fasta file. Returns the
    path to the newly created kallisto index as a string.
    '''
    output_file = os.path.join(output_dir, kalliso_index_name)
    cmd = [kallisto_executable, 'index', '-i', output_file, trans_file,
           '--make-unique']
    subprocess.call(cmd)
    return output_file


