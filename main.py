import os
import sys

from bowtie import build_bowtie_index, search_BTI
from data import (convert_to_fastq, download_accession, get_files_from_parent,
                  get_paired_end_paths_as_lists, get_paired_paths_outer_dir,
                  run_wget)
from get_args import get_args
from kallisto import make_kalisto_index, run_kallisto
from qBLAST import run_and_write_blast
from sleuth import make_sleuth_table, run_sleuth
from spades import (assemble_with_spades, concat_contigs, make_big_fasta,
                    write_assembly_stats)


def main():

    args = get_args()
    log = open(os.path.join(args.o, 'miniproject.log'), 'w')
    log.write('Comp Bio Mini Project\n')

    if not args.q:  # dont need paths to SRA files if already have fastq
        if args.i:  # already have the input files
            args.i = [os.path.join(args.i, f) for f in os.listdir(args.i)]
        else:
            args.i = run_wget(args.o)  # right now would return with file name

    # now have abs paths to all files at this point
    # files in fastq_format
    print('Downloading accessions')
    trans_file_path = download_accession(args.o, log)
    # get CDS and put into fasta file
    print('Making kalliso index')
    if not args.k:  # kallisto index already created can skip making it
        args.k = make_kalisto_index(trans_file_path, args.o)

    if not args.q:
          # no fastq files given so convert SRA files in input dir
        args.q = get_paired_end_paths_as_lists(
            convert_to_fastq(args.i, args.o))
    else:
        args.q = get_paired_paths_outer_dir(args.q)

    if not args.r:
        args.r = run_kallisto(args.k, args.q, args.o)
    else:
        args.r = get_files_from_parent(args.r)

    sleuth_table = make_sleuth_table(args.r, args.o)
    print('Running Sleuth')
    sleuth_results = run_sleuth(sleuth_table, args.o, log)
    print('Downloading Complete Genome')
    complete_genome = download_accession(args.o, log, dtype='genome')
    print('Running Bowtie2')
    if not args.b:
        args.b = build_bowtie_index(complete_genome, args.o)

    if not args.s:
        args.s = search_BTI(args.b, args.q, args.o, log, threads=args.t)
    else:
        args.s = get_files_from_parent(args.s)  # these are sam files

    if not args.f:  # if given would still need to convert them
        args.f = make_big_fasta(args.s, args.o)

    if not args.a:
        args.a = assemble_with_spades(args.f, args.o, args.t)

    print('Writing Assembly Stats to log file')
    # pass in contigs file path and write stats to log
    write_assembly_stats(args.a, log)
    cat = concat_contigs(args.a)
    print('Running BLAST, query length = {}'.format(len(cat)))
    print(args.local)
    run_and_write_blast(cat, args.o, log, local=args.local)
    log.close()
    print('==============================\nRun Complete')


if __name__ == '__main__':
    main()
