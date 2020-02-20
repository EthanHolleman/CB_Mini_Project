import subprocess
import os

from Bio import Entrez
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq


def run_kallisto(k_index, query_file_paths, kallisto_executable='kallisto'):
    pass
    '''
    Given the path to a kallisto index and a list of query files
    takes the query filepaths uses subprocess to run analysis using 
    kallisto program. 
    '''


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


def download_accession(output_dir, entrez_email='eholleman@luc.edu'):
    '''
    Use Biopython Entrez and SeqIO to first pull the transcriptome accession
    from genbank(nucleotide database) and then write the records in the
    SeqIO object to a fasta file which then can be passed into the 
    make_kalisto_index function. Returns the absolute path to the 
    accession fasta file as a string and the number of coding sequences 
    as an int both in a tuple.
    '''
    ACC = 'EF999921.1'
    output_file = os.path.join(output_dir, ACC)
    Entrez.email = entrez_email
    handle = Entrez.efetch(db="nucleotide", id=ACC,
                           rettype="gb", retmode="text")
    #  fetch the accession
    record = SeqIO.read(handle, 'genbank')
    # parse and a genbank file
    record.features = [f for f in record.features if f.type == "CDS"]
    num_cds = len(record.features)
    # make sure features only contain the CDS sequences
    seq_recs = [r.extract(record) for r in record.features]
    # convert the SeqFeatures objects to reqrecords so can write using SeqIO
    SeqIO.write(seq_recs, output_file, 'fasta')  # write as Fasta file

    return output_file, num_cds
