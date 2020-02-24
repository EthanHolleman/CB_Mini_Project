import subprocess
import os

from Bio.Seq import Seq
from Bio import SeqIO

from data import get_files_from_parent


def assemble_with_spades(input_file, output_dir, threads=4,
                         spades_executable='spades'):
    '''
    Reads in a list of fasta files and assembles using spades. Returns the
    output dir path. 
    '''

    cmd = [spades_executable, '--12', input_file, '-k', '55,77,99,127', '-t', threads,
           '--only-assembler', '-o', output_dir]
    cmd = [str(c) for c in cmd]  # convert everything to string
    subprocess.call(cmd)

    return output_dir


def make_big_fasta(fasta_files, output_dir, big_fasta_name='big_fasta.fa'):
    '''
    Given a list of fasta files created from the converting the sam file
    output of a bowtie alignment uses cat to create one large fasta
    file to pass onto spades for assembly.
    '''
    big_fasta_name = os.path.join(output_dir, big_fasta_name)
    cmd = ' '.join(['cat'] + fasta_files + ['>', big_fasta_name])
    subprocess.call(cmd, shell=True)

    return big_fasta_name


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
