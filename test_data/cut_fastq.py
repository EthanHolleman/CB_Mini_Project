import os
from subprocess import check_output
import math

SRA_to_FASTQ = '/home/ethan/Documents/kallisto_test/SRA_to_FASTQ'
PROP = 0.25

def round_up_to_even(f):
    return math.ceil(f / 2.) * 2

def get_fastx_length(fastx, t='q'):
    cmd = ['wc', '-l', fastx]
    call = check_output(cmd)
    lines = call.decode('utf-8').split(' ')[0]
    if t == 'q':
        return int(lines) / 4   # if fastq 4 lines per entry
    elif t == 'a':
        return int(lines) / 2  # fasta 2 lines per entry


sub_dirs = [os.path.join(SRA_to_FASTQ, d) for d in os.listdir(SRA_to_FASTQ)]

for sub_dir in sub_dirs:
    fastq_files = [os.path.join(sub_dir, f) for f in os.listdir(sub_dir)]
    for fq in fastq_files:
        test_name = fq[:-6] + '_test.fastq'
        l = get_fastx_length(fq)
        nl = round_up_to_even(get_fastx_length(fq) * 0.05)
        os.system('head {} -n {} > {}'.format(fq, nl, test_name))
    