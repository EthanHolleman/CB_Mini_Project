import subprocess
import os

URLS = ['https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1',
        'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1',
        'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1',
        'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1']


def run_wget(output_dir, urls=URLS):
    paths = []
    for url in urls:
        p = os.path.join(output_dir, os.path.basename(url))
        cmd = ['wget', url, '-O', p]
        c = subprocess.call(cmd)
        paths.append(p)

    return paths


def convert_to_fastq(SRA_paths):
    for SRA in SRA_paths:
        cmd = ['fastq-dump', SRA]
        subprocess.call(cmd)
        yield SRA_paths + '.fastq'
        #  yield the fastq file path
