import subprocess
import os

URLS = {'Donor_1_first': 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1',
'Donor_1_after': 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1',
'Donor_3 first': 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1',
'Donor 3 after': 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1'}


def run_wget(urls, output_dir):
    for url in urls:
        cmd = ['wget', url, output_dir]
        c = subprocess.call(cmd)

    return [os.path.join(output_dir, f) for f in os.listdir(output_dir)]

def convert_to_fastq(SRA_paths):
    for SRA in SRA_paths:
        cmd = ['fastq-dump', SRA]
        subprocess.call(cmd)
        yield SRA_paths + '.fastq'
        #  yield the fastq file path
