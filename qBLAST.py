import csv
import os

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import Seq

# run all blasting operations from this file


def run_blast(seq_object):
    '''
    Given a string of a seq object blasts against nucleotide database
    for Herp only. Returns a blast record type object object
    that can be parsed using NCBIXML.read()
    '''
    return NCBIWWW.qblast('blastn', 'nr', str(seq_object), entrez_query='Herpesviridae[ORGN]')
    # use an entrez query to set limit the search range


def get_top_ten_results(xml_file):
    '''
    Given XML file returns the top ten results of the blast query.
    Assumes that results are list in accesnding order according to
    their evalues. Returns a list for each result that can then
    be written using csv writter to log file.
    '''
    handle, top_hits = open(xml_file), []
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


def write_top_hits(top_hits, out_dir, file_name='top_10.tsv'):  # probably want to replace this with log file name
    outfile = os.path.join(out_dir, file_name)
    with open(outfile, 'w') as of:
        writer = csv.writer(of, delimiter='\t')
        for row in top_hits:
            writer.writerow(row)
    return outfile

