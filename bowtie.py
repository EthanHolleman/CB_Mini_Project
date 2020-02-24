import subprocess
from subprocess import check_output
import os

from Bio import SeqIO
from Bio.Seq import Seq
from data import if_not_dir_make

from sleuth import default_conditions, test_for_condition

DONOR_DICT = {'SRR5660030.1': 'Donor 1', 'SRR5660033.1': 'Donor 1',
              'SRR5660044.1': 'Donor 3', 'SRR5660045.1': 'Donor 3'}


def build_bowtie_index(input_file, output_dir, index_name='MP_BTI'):
    '''
    Given an index name, input fasta file and an output path creates
    a new bowtie2 index and returns the path to that index for querying.
    '''
    BT_path = if_not_dir_make(output_dir, index_name)
    index_path = os.path.join(BT_path, index_name)
          
    print('Building Bowtie Index')
    cmd = ['bowtie2-build', input_file, index_path]
    subprocess.call(cmd)

    return index_path

def search_BTI(BTI, query_files, output_dir, log, threads=4):
    '''
    Given a bowtie2 index and a list of list of paired end reads
    where each sublist has the paired reads. Alligns the reads
    against the bowtie index and then converts to fastq format.

    a query file alligns sequnces in the query file
    (fasta format) to the index and writes the allignment SAM file to the
    output dir.
    '''
    cond_dict = default_conditions()  # use to test for which condition
    output_dir = if_not_dir_make(output_dir, 'bowtie_results')
    
    
    output_files = []
    for read_a, read_b in query_files:
        sam_name = os.path.basename(read_a) + '.sam'
        output_file = os.path.join(output_dir, sam_name)
        
        sample, date = test_for_condition(cond_dict, read_a)
        donor = DONOR_DICT[sample]
        
        print('Running new Bowtie Search')
        cmd = ['bowtie2', '-x', BTI, '-1', read_a, '-2',
               read_b, '-S', output_file, '--threads', '4', '--no-discordant', 
               '--no-unal', '--no-mixed']
        print(' '.join(cmd))
        subprocess.call(cmd)
        fasta_file = convert_sam_to_fasta(output_file)
        output_files.append(fasta_file)
        log_string = make_read_comparison_string(read_a, fasta_file, donor, date)
        log.write(log_string + '\n')
        print(log_string)

    return output_files


def convert_sam_to_fasta(sam_file):
    '''
    Uses sam tools to convert a sam file from bowtie to a fasta
    formated file. Returns the path to the new fasta file.
    '''
    fasta_name = sam_file + '.fasta'
    cmd = ' '.join(['samtools', 'fasta', sam_file, '>', fasta_name])
    os.system(cmd)

    return fasta_name


def make_read_comparison_string(intial, postBT, donor, date):
    '''
    Given a fastq file before alignment using bowtie counts the
    number of reads in that file and in the fasta file created
    from its sam file output. Assumes that the reads in the fasta
    file are paired end.
    ''' 
    # divide by 2 again becuase interlaced paired ends
    i, p = get_fastx_length(intial), get_fastx_length(postBT, 'a') / 2
    return '{} ({}) had {} reads pairs before Bowtie2 filtering and \
            {} read pairs after.'.format(donor, date, i, p.strip())


def get_fastx_length(fastx, t='q'):
    cmd = ['wc', '-l', fastx]
    call = check_output(cmd)
    lines = call.decode('utf-8').split(' ')[0]
    if t == 'q':
        return int(lines) / 4   # if fastq 4 lines per entry
    elif t == 'a':
        return int(lines) / 2  # fasta 2 lines per entry
