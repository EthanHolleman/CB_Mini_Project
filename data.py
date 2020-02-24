import subprocess
import os

from Bio import Entrez
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq

URLS = ['https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1',
        'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1',
        'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1',
        'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1']

# store urls in a constant 

def if_not_dir_make(parent_dir, dir_name):
    full_dir = os.path.join(parent_dir, dir_name)
    if not os.path.isdir(full_dir):
        os.mkdir(full_dir)
    return full_dir


def run_wget(output_dir, urls=URLS):
    '''
    Given a list of valid urls (file locations) uses subprocess to run wget
    to download all files in that list of urls. Writes downloaded files to
    the given output_dir and names files based on the basename of the url
    (everything after the last /). Returns a list of absoulte paths to the
    downloaded files as a list. 
    '''
    paths = []
    for url in urls:
        p = os.path.join(output_dir, os.path.basename(url))
        cmd = ['wget', url, '-O', p]
        subprocess.call(cmd)
        paths.append(p)

    return paths


def convert_to_fastq(SRA_paths, output_dir):
    '''
    Given a list of paths so SRA_paths that should have been downloaded using
    the run_wget function uses fastq-dump command line program through 
    subprocess to convert all SRA files to fastq files. Yeilds paths to the
    fastq files as it runs. Appends ''
    for SRA in SRA_paths appends .fastq to the filenames.
    '''
    print('Converting files to fastq format')
    
    output_dir = if_not_dir_make(output_dir, 'SRA_to_fastq')

    # need to get paths to paired end files
    fastq_paths = []
    for SRA in SRA_paths:
        fastq_name = os.path.join(output_dir, os.path.basename(SRA) + '.fastq')
        cmd = ['fastq-dump', SRA, '-O', fastq_name, '--split-files']
        subprocess.call(cmd)
        # returns a directory with the paired end files
        fastq_paths.append(fastq_name)
    return fastq_paths


def get_paired_end_paths_as_lists(fastq_dirs):
    '''
    Given the output from convert_to_fastq (directories holding the newly
    converted paired end fastq files) returns as list of lists the paths to
    the individual fastq files so they can be passed into kallisto function.
    '''
    paired_reads = []
    for fq_dir in fastq_dirs:
        paired_reads.append([os.path.join(fq_dir, f)
                             for f in os.listdir(fq_dir)])
    return paired_reads


def get_paired_paths_outer_dir(fastq_dir):
    '''
    Given a parent dir that contains seperate dirs for each set of
    paired end reads returns a list of lists where each sublist contains the
    two paired end fastq files.
    '''
    paired_reads = []
    fastq_dirs = [os.path.join(fastq_dir, q) for q in os.listdir(
        fastq_dir)]  # get all fastq containing dirs
    for fq_dir in fastq_dirs:
        paired_reads.append([os.path.join(fq_dir, f)
                             for f in os.listdir(fq_dir)])
    return paired_reads


def download_accession(output_dir, entrez_email='eholleman@luc.edu',
                       dtype='cdna'):
    '''
    Use Biopython Entrez and SeqIO to first pull the transcriptome accession
    from genbank(nucleotide database) and then write the records in the
    SeqIO object to a fasta file which then can be passed into the 
    make_kalisto_index function. Returns the absolute path to the 
    accession fasta file as a string and the number of coding sequences 
    as an int both in a tuple. 

    Download accession also works in another mode. By setting the dtype arg to 
    'genome' it will download and return path to the complete genome of the 
    accession. Used for handing off to bowtie index maker. 
    '''
    ACC = 'EF999921.1'
    output_file = os.path.join(output_dir, ACC)
    Entrez.email = entrez_email
    handle = Entrez.efetch(db="nucleotide", id=ACC,
                           rettype="gb", retmode="text")
    #  fetch the accession
    record = SeqIO.read(handle, 'genbank')

    if dtype == 'cdna':  # for just coding sequences
        record.features = [f for f in record.features if f.type == "CDS"]
        num_cds = len(record.features)  # store for log
        seq_recs = [r.extract(record) for r in record.features]
        SeqIO.write(seq_recs, output_file + '_cdna', 'fasta')

        return output_file, num_cds
    elif dtype == 'genome':  # for whole genome 
        SeqIO.write(record, output_file + '_genome', 'fasta')

        return output_file
    
    
def get_files_from_parent(parent_dir):
    '''
    General helper function to return the absolute paths of all files in
    a parent directory. Used for passing in directory of files instead of
    the individual file paths.
    '''
    return [os.path.join(parent_dir, d) for d in os.listdir(parent_dir)]

  