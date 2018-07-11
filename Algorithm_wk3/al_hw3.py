from editDistance import *
from Assembly import *

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def readFASTQ(filename):
	sequences = []
	qualities = []
	with open (filename) as f:
		while True:
			# skip name line
			f.readline()
			# store sequence of the read
			seq = f.readline().rstrip()
			# skip spaceholder line
			f.readline()
			# store qulity score of the read
			qual = f.readline().rstrip()
			# if we reach the end of the line,break and finish
			if len(seq) == 0:
				break
			sequences.append(seq)
			qualities.append(qual)
		return sequences, qualities


# test case 2
reads = ['ABCDEFG', 'EFGHIJ', 'HIJABC']
olaps,outgoing = overlap_all_pairs(reads, 3)
print len(olaps)

reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
olaps,outgoing = overlap_all_pairs(reads, 4)
print len(olaps)

reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
olaps,outgoing = overlap_all_pairs(reads, 5)
print len(olaps)

# homework
t = readGenome('chr1.GRCh38.excerpt.fasta')

p = 'GCTGATCGATCGTACG'
editDist = editDistance_modified(p, t)
print 'Q1: the edit distance for the best match is: %d' % editDist

p = 'GATTTACCAGATTGAG'
editDist = editDistance_modified(p, t)
print 'Q2: the edit distance for the best match is: %d' % editDist


reads, quals = readFASTQ('ERR266411_1.for_asm.fastq')

olaps,outgoing = overlap_all_pairs(reads,30)
print 'Q3: # of overlaps are: %d' % len(olaps)
print 'Q4: # of outgoing nodes are: %d' %len(outgoing)
