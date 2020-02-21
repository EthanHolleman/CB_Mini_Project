from Bio.Blast import NCBIWWW

# run all blasting operations from this file

def run_blast(seq_object):
    return NCBIWWW.qblast('blastn', 'nr', seq_object, entrez_query='Herpesviridae[ORGN]')
	# use an entrez query to set limit the search range
 