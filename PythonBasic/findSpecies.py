from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
fasta_string = """TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTAC
AATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCAC
CTACGGTAGAG"""
result_handle = NCBIWWW.qblast("blastn","nt",fasta_string)

blast_record = NCBIXML.read(result_handle)

E_VALUE_THRESH = 0.01
for alignment in blast_record.alignments:
	for hsp in alignment.hsps:
		if hsp.expect<E_VALUE_THRESH:
			print "***Alignment***"
			print "sequence:", alignment.title
			print "length", alignment.length
			print "e evalue", hsp.expect
			print hsp.query
			print hsp.match
			print hsp.sbjct
