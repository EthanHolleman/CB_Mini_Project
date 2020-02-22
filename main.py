import os
from get_args import get_args
from data import run_wget, convert_to_fastq, get_paired_end_paths_as_lists, download_accession
from kallisto import make_kalisto_index, run_kallisto
from sleuth import make_sleuth_table

def main():

    args = get_args()
    if args.i:  # already have the input files
        args.i = [os.path.join(args.i, f) for f in os.listdir(args.i)]
    else:
        args.i = run_wget(args.o)  # right now would return with file name

    # now have abs paths to all files at this point
    # files in fastq_format
    trans_file_path, num_cds = download_accession(args.o)
    # get CDS and put into fasta file
    if not args.k:  # kallisto index already created can skip making it
        args.k = make_kalisto_index(trans_file_path, args.o)
        
    if not args.q:
        args.q = convert_to_fastq(args.i)
    
    paired_fastqs = get_paired_end_paths_as_lists(args.q)
    
    kallisto_dirs = run_kallisto(args.k, paired_fastqs, args.o)
    
    sleuth_table = make_sleuth_table(kallisto_dirs, args.o)
    
    
    
    
    #print(args.q)

    #args.q = get_paired_end_paths_as_lists(args.i)
    
    #run_kallisto(args.k, args.q, args.o)
    # run kallisto program on all fastq files need to look into passing
    # args into this for keeping track of what files go with what
    # maybe include some type of config file to do this?

    # kallisto index is made and files are ready to go at this point
    # issue if already have index donesnt download the
    # need to write below to log file
    # The HCMV genome (EF99921) has #CDS


if __name__ == '__main__':
    main()
