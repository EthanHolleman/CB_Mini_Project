import subprocess
import os

URLS = ['https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1',
        'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1',
        'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1',
        'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1']


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


def convert_to_fastq(SRA_paths):
    '''
    Given a list of paths so SRA_paths that should have been downloaded using
    the run_wget function uses fastq-dump command line program through 
    subprocess to convert all SRA files to fastq files. Yeilds paths to the
    fastq files as it runs. Appends ''
    for SRA in SRA_paths appends .fastq to the filenames.
    '''
    print('Converting files to fastq format')
    
    # need to get paths to paired end files
    fastq_paths = []
    for SRA in SRA_paths:
        fastq_name = SRA + '_fastq'
        cmd = ['fastq-dump', SRA, '-O', fastq_name, '--split-files']
        subprocess.call(cmd)
        fastq_paths.append(fastq_name)  # returns a directory with the paired end files
    return fastq_paths

def get_paired_end_paths_as_lists(fastq_dirs):
    '''
    Given the output from convert_to_fastq (directories holding the newly
    converted paired end fastq files) returns as list of lists the paths to
    the individual fastq files so they can be passed into kallisto function.
    '''
    paired_reads = []
    for fq_dir in fastq_dirs:
        paired_reads.append([os.path.join(fq_dir, f) for f in os.listdir(fq_dir)])
    return paired_reads
