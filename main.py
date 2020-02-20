import os
from get_args import get_args
from data import run_wget, convert_to_fastq
from kallisto import download_accession, make_kalisto_index, run_kallisto \
    convert_to_fastq



def main():

    args = get_args()
    if args.i:  # already have the input files
        SRA_paths = [os.path.join(args.i, f) for f in os.listdir(args.i)]
    else:
        SRA_paths = run_wget(args.o)  # right now would return with file name

    # now have abs paths to all files at this point
    # files in fastq_format
    fastq_paths = convert_to_fastq(SRA_paths)
    trans_file_path, num_cds = download_accession(args.o)
    # get CDS and put into fasta file
    if not args.k:  # kallisto index already created can skip making it
        args.k = make_kalisto_index(trans_file_path, args.o)
    
    if not args.q:  # SRA files have not been converted to fastq
        args.q = convert_to_fastq(SRA_paths)
    
    
        
        
    
    run_kallisto(args.k, fastq_paths)
    # run kallisto program on all fastq files need to look into passing
    # args into this for keeping track of what files go with what
    # maybe include some type of config file to do this?
    

    # kallisto index is made and files are ready to go at this point
    # issue if already have index donesnt download the
    # need to write below to log file
    # The HCMV genome (EF99921) has #CDS

if __name__ == '__main__':
    main()