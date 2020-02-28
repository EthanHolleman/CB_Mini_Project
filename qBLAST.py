import csv
import os
import subprocess

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO

from data import if_not_dir_make

# run all blasting operations from this file


def run_and_write_blast(query, out_dir, log, local=0):
    '''
    Function that wraps up run_blast, get_top_results and
    write_top_results into one function. Takes in a query string
    an output dir and the log file. Runs a blast search using biopython
    and writes the results to a new subdir called BLAST_results in an 
    xml file. Then takes the top ten results from that xml file and
    writes them to log file along with a header in tsv format.
    '''
    if local == 1 or local == '1':
        bdb = make_local_BDB()
        xml_file = local_blast(query, bdb, out_dir)
    else:
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
    xml_path = os.path.join(blast_path, dir_name + '.xml')
    xml = NCBIWWW.qblast('blastn', 'nr', str(
        seq_object), entrez_query='Herpesviridae[ORGN]', hitlist_size=10, expect=1e-200, megablast=True, alignments=10)

    with open(xml_path, "w") as out_handle:
        print('writing results')
        out_handle.write(xml.read())

    return xml_path
    # use an entrez query to set limit the search range


def get_top_ten_results(xml_file_path):
    '''
    Given XML file returns the top ten results of the blast query.
    Assumes that results are list in accesnding order according to
    their evalues. Returns a list for each result that can then
    be written using csv writter to log file.
    '''
    with open(xml_file_path) as handle:
        top_hits, i = [], 0
        blast_record = NCBIXML.read(handle)
        while i < 10:
            alignment = blast_record.alignments[i]
            hsp = alignment.hsps[0]
            values = [alignment.title, alignment.length,
                      len(alignment.hsps), hsp.identities,
                      hsp.gaps, hsp.bits, hsp.expect]
            for j, v in enumerate(values):
                if v == None:
                    values[j] = 0
            top_hits.append(values)
            i += 1
        return top_hits


def write_top_hits(top_hits, log):  # probably want to replace this with log file name
    '''
    Writes top hits from list of lists to the given log file in tsv format.
    '''
    HEADER = [['seq_title', 'align_len', 'number_HSPs', 'topHSP_ident',
               'topHSP_gaps', 'topHSP_bits', 'topHSP_expect']]
    writer = csv.writer(log, delimiter='\t')
    writer.writerows(HEADER + top_hits)


def local_blast(seq, local_BDB, output_dir, dir_name='BLAST_results'):
    blast_path = if_not_dir_make(output_dir, dir_name)
    xml_path = os.path.join(blast_path, dir_name + '.xml')
    query_path = os.path.join(blast_path, 'BLAST_query.fasta')
    print('writing query')
    print(query_path)
    with open(query_path, 'w') as qp:
        qp.write('>Query Sequence\n')
        qp.write(seq)
    cmd = ['blastn', '-db', local_BDB, '-query',
           query_path, '-outfmt', '5', '-out', xml_path]
    subprocess.call(cmd)

    return xml_path


def make_local_BDB(xz_path='./test_data/Herp_BDB.tar.xz', xz_dir='./test_data', BDB_name='Herp_BDB'):
    cmd = ['tar', 'xf', xz_path, '-C', xz_dir]
    subprocess.call(cmd)
    db_name = os.path.join(xz_dir, BDB_name)
    print(db_name)

    return db_name
