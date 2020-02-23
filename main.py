import os
from get_args import get_args
from data import run_wget, convert_to_fastq, get_paired_end_paths_as_lists, download_accession
from data import get_paired_paths_outer_dir
from kallisto import make_kalisto_index, run_kallisto
from kallisto import get_kallisto_dir_paths
from sleuth import make_sleuth_table, run_sleuth
from bowtie import build_bowtie_index
from bowtie import search_BTI

def main():

    args = get_args()
    
    if not args.q:  # dont need paths to SRA files if already have fastq
        if args.i:  # already have the input files
            args.i = [os.path.join(args.i, f) for f in os.listdir(args.i)]
        else:
            args.i = run_wget(args.o)  # right now would return with file name

    # now have abs paths to all files at this point
    # files in fastq_format
    trans_file_path, num_cds = download_accession(args.o)
    print('downloaded acc')
    # get CDS and put into fasta file
    print('making index')
    if not args.k:  # kallisto index already created can skip making it
        args.k = make_kalisto_index(trans_file_path, args.o)
    

    if not args.q:
          # no fastq files given so convert SRA files in input dir
        args.q = get_paired_end_paths_as_lists(convert_to_fastq(args.i) )
    else:
        args.q = get_paired_paths_outer_dir(args.q)
    
    if not args.r:
        args.r = run_kallisto(args.k, args.q, args.o)
    else:
        args.r = get_kallisto_dir_paths(args.r)
    
    sleuth_table = make_sleuth_table(args.r, args.o)
    print('Running Sleuth')
    sleuth_results = run_sleuth(sleuth_table, args.o)
    
    complete_genome = download_accession(args.o, dtype='genome')
    
    if not args.b:
        args.b = build_bowtie_index(complete_genome, args.o)
    
    sam_files = search_BTI(args.b, args.q, args.o)
    
        
    
    
    
    
    # run the sleuth shit and write that to a place to read into by next functions
    
    
    
    
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
