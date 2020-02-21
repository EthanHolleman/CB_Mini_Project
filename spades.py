import subprocess

from Bio.Seq import Seq
from Bio import SeqIO


def assemble_with_spades(input_files, output_dir, threads=4,
                         spades_executable='spades'):
    '''
    Reads in a list of fasta files and assembles using spades. Returns the
    output dir path. 
    '''
    input_list = ['-{} {}'.format(i+1, f) for i, f in enumerate(input_files)]
    cmd = [spades_executable, '-k', '55,77,99,127', '-t', threads,
           '--only-assembler', '-o', output_dir] + input_list
    subprocess.call(cmd)

    return output_dir

# outputs a contigs and scaffolds fasta files in the output dir


def read_fasta_with_bio(fasta_file):
    return SeqIO.parse(fasta_file, 'fasta')


def count_contigs(fasta_file, min_length=1000):
    '''
    Given a fasta file returns a count of the number of sequences that have
    a length greater than min_length. Defualt min_length = 1000.
    '''
    records = read_fasta_with_bio(fasta_file)

    return len([r for r in records if len(r.seq) > min_length])


def count_length_assembly(fasta_file):
    '''
    Given the path to a fasta file returns the total number of base pairs
    of all the sequences in the file.
    '''
    records = read_fasta_with_bio(fasta_file)
    return sum([len(r.seq) for r in records])


def concat_contigs(fasta_file, min_length):
    '''
    Given a fasta file concatinated all the sequences with length greater than
    the min_length arg, defualt 1000. Adds 50 Ns inbetween each concatenated
    sequence.
    '''
    records = read_fasta_with_bio(fasta_file)
    cat_seq = []
    for r in records:
        if len(r.seq) > 1000:
            cat_seq.append(str(r.seq))
            cat_seq.append('N'*50)
    return ''.join(cat_seq)
