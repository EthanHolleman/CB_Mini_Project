import subprocess
import os


def build_bowtie_index(input_file, output_dir, index_name='MP_BTI'):
    '''
    Given an index name, input fasta file and an output path creates 
    a new bowtie2 index and returns the path to that index for querying.
    '''
    BTI = os.path.join(output_dir, index_name)
    print('Building Bowtie Index')
    cmd = ['bowtie2-build', input_file, BTI]
    subprocess.call(cmd)

    return BTI


def search_BTI(BTI, query, output_dir, sam_name='EF999921.sam'):
    '''
    Given a bowtie2 index and a query file alligns sequnces in the query file
    (fasta format) to the index and writes the allignment SAM file to the
    output dir. 
    '''
    output_file = os.path.join(output_dir, sam_name)
    cmd = ['bowtie2', BTI, '-f', query, '-S', output_file]
    subprocess.call(cmd)
    
    return output_file