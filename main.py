import os
from get_args import get_args
from data import run_wget, convert_to_fastq
from kallisto import download_accession, make_kalisto_index

def main():

    args = get_args()
    if args.i:  # already have the input files
        SRA_paths = [os.path.join(args.i, f) for f in os.listdir(args.i)]
    else:
        SRA_paths = run_wget()  # right now would return with file name

    # now have abs paths to all files at this point
    # files in fastq_format
    fastq_paths = convert_to_fastq(SRA_paths)
    trans_file_path, num_cds = download_accession(args.o)
    # get CDS and put into fasta file
    if not args.k:  # kallisto index already created can skip making it
        args.k = make_kalisto_index(trans_file_path, args.o)

    # kallisto index is made and files are ready to go at this point
    # issue if already have index donesnt download the
    # need to write below to log file
    # The HCMV genome (EF99921) has #CDS

    
