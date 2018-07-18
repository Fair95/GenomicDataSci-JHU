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


# test case 1
strings = ['ABC', 'BCA', 'CAB']
print scs_list(strings)
strings = ['GAT', 'TAG', 'TCG', 'TGC', 'AAT', 'ATA']
# Returns list of all superstrings that are tied for shorest
print scs_list(strings)

ss = ['CCT','CTT','TGC','TGG','GAT','ATT']
print 'Q1: length of the scs is %d' % len(scs(ss))
print 'Q2: # of different scs is %d' % len(scs_list(ss))

reads,quals = readFASTQ('ads1_week4_reads.fq')
assemble = greedy_scs_optimized(reads,10)
print assemble.count('A')
print assemble.count('T')