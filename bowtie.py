import subprocess
from subprocess import check_output
import os

from Bio import SeqIO
from Bio.Seq import Seq


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

# need to do a paired end search


def search_BTI(BTI, query_files, output_dir, sam_name='EF999921.sam', threads=4):
    '''
    Given a bowtie2 index and a list of list of paired end reads
    where each sublist has the paired reads. Alligns the reads
    against the bowtie index and then converts to fastq format.
    
    a query file alligns sequnces in the query file
    (fasta format) to the index and writes the allignment SAM file to the
    output dir.
    '''

    output_files = []
    for read_a, read_b in query_files:
        sam_name = os.path.basename(read_a) + '.sam'
        query_file_len = get_fastq_length(read_a)
        output_file = os.path.join(output_dir, sam_name)
        output_files.append(output_file)
        print('Running new Bowtie Search')
        cmd = ['bowtie2', '-x', BTI, '-1', read_a, '-2', read_b, '-S', output_file, '--threads', '4']
        subprocess.call(cmd)

    return output_files


r = '/media/ethan/KINGSTON/kallisto_test/SRA_to_fastq/SRR5660044.1_fastq/SRR5660044.1_2.fastq'


def get_fastq_length(fastq):
    cmd = ['wc', '-l', fastq]
    call = check_output(cmd)
    return call.decode('utf-8').split(' ')[0]

