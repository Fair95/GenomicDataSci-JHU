def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def strandAwareNaive(p,t):
    if p == reverseComplement(p):
        occurrences = naive(p,t)
        return occurrences
    else:
        occurrences = naive(p,t)
        occurrences.extend(naive(reverseComplement(p),t))
        return occurrences


def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        mismatch = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mismatch += 1
            if mismatch >2:
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

lambda_genome = readGenome('lambda_virus.fa')
occurrences = strandAwareNaive('AGGT', lambda_genome)
print('Q1: # occurrences: %d' % len(occurrences))

occurrences = strandAwareNaive('TTAA', lambda_genome)
print('Q2: # occurrences: %d' % len(occurrences))

occurrences = strandAwareNaive('ACTAAGT', lambda_genome)
print('Q3: offset of leftmost occurrence: %d' % min(occurrences))

occurrences = strandAwareNaive('AGTCGA', lambda_genome)
print('Q4: offset of leftmost occurrence: %d' % min(occurrences))

occurrences = naive_2mm('TTCAAGCC', lambda_genome)
print('Q5: # occurrences: %d' % len(occurrences))

occurrences = naive_2mm('AGGAGGTT', lambda_genome)
print('Q6: offset of leftmost occurrence: %d' % min(occurrences))

test_genome,qualities = readFastq("ERR037900_1.first1000.fastq")
from algorithm_wk1 import Phred33ToQ
quals = ''.join(qualities)
for i in range(len(quals)):
    if Phred33ToQ(quals[i]) == 2:
        print "Q7: offset of leftmost bad quality circle: %d" %i
        break
