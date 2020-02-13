import subprocess
import os

URLS = ['https://www.ncbi.nlm.nih.gov/sra/SRX2896360',
'https://www.ncbi.nlm.nih.gov/sra/SRX2896363',
'https://www.ncbi.nlm.nih.gov/sra/SRX2896374',
'https://www.ncbi.nlm.nih.gov/sra/SRX2896375']


def run_wget(urls, output_dir):
    for url in urls:
        cmd = ['wget', url, output_dir]
        c = subprocess.call(cmd)

    return [os.path.join(output_dir, f) for f in os.listdir(output_dir)]
