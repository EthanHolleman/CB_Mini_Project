import csv
import os

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import Seq

from data import if_not_dir_make

# run all blasting operations from this file


def run_and_write_blast(query, out_dir, log):
    xml_file = run_blast(query, out_dir)
    top_results = get_top_ten_results(xml_file)
    write_top_hits(top_results, log)


def run_blast(seq_object, output_dir, dir_name='BLAST_results'):
    '''
    Given a string of a seq object blasts against nucleotide database
    for Herp only. Returns a blast record type object object
    that can be parsed using NCBIXML.read()
    '''
    blast_path = if_not_dir_make(output_dir, dir_name)
    xml = NCBIWWW.qblast('blastn', 'nr', str(seq_object), entrez_query='Herpesviridae[ORGN]')
    
    with open(blast_path, "w") as out_handle:
        out_handle.write(xml.read())
    
    return xml
    # use an entrez query to set limit the search range


def get_top_ten_results(xml_file_path):
    '''
    Given XML file returns the top ten results of the blast query.
    Assumes that results are list in accesnding order according to
    their evalues. Returns a list for each result that can then
    be written using csv writter to log file.
    '''
    with open(xml_file_path) as handle:
        top_hits = []
        blast_record = NCBIXML.read(handle)
        for i, alignment in enumerate(blast_record.alignments):
            if i == 10:
                break
            for hsp in alignment.hsps:
                if i == 0:
                    top_hits.append([alignment.title, alignment.length,
                                    hsp.num_alignments, hsp.identities,
                                    hsp.gaps, hsp.bits, hsp.expect])
                else:
                    top_hits.append([alignment.title, alignment.length])
        return top_hits


def write_top_hits(top_hits, log):  # probably want to replace this with log file name
    writer = csv.writer(log, delimiter='\t')
    for row in top_hits:
        writer.writerow(row)

#from spades import concat_contigs
#f = '/home/ethan/Documents/spades_contigs/contigs.fasta'
#string = concat_contigs(f)
#print(string)
#print('Running Blast')
#xml = run_blast(string, './')

#print('Finished Blast')
#print(get_top_ten_results(xml))
