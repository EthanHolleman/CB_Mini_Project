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


def search_BTI(BTI, query_files, output_dir, threads=4):
    '''
    Given a bowtie2 index and a list of list of paired end reads
    where each sublist has the paired reads. Alligns the reads
    against the bowtie index. Results are returned as a sam file to the
    location given by the output_dir variable.
    '''
    output_files = []
    for read_a, read_b in query_files:
        sam_name = os.path.basename(read_a) + '.sam'
        query_file_len = get_fastq_length(
            read_a)  # store len for writing to log
        output_file = os.path.join(output_dir, sam_name)
        output_files.append(output_file)  # store path to output
        print('Running new Bowtie Search')
        cmd = ['bowtie2', '-x', BTI, '-1', read_a, '-2',
               read_b, '-S', output_file, '--threads', '4']
        subprocess.call(cmd)

    return output_files


def get_fastq_length(fastq):
    '''
    Helper command that uses wc -l to get the number of entries in a fastq
    file. Each entry should have four lines so divides the results of the
    wc call by 4 to get the number of entries. Paired end read files should
    have the same number of entries so only needs to be called on one of the
    mates.
    '''
    cmd = ['wc', '-l', fastq]
    call = check_output(cmd)
    return int(call.decode('utf-8').split(' ')[0]) / 4
