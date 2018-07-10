from bm_preproc import *
from pigeonhole_principle import *
import wget
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def naive_with_counts(p, t):
    occurrences = []
    num_character_comparisons = 0
    num_alignments = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        num_alignments += 1
        for j in range(len(p)):  # loop over characters
            num_character_comparisons += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, num_alignments, num_character_comparisons


         
# testing case 1
p = 'word'
t = 'there would have been a time for such a word'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print(occurrences, num_alignments, num_character_comparisons)

p = 'needle'
t = 'needle need noodle needle'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print(occurrences, num_alignments, num_character_comparisons)

p = 'word'
t = 'there would have been a time for such a word'
lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print(occurrences, num_alignments, num_character_comparisons)

p = 'needle'
t = 'needle need noodle needle'
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print(occurrences, num_alignments, num_character_comparisons)

## testing case 2
t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
p = 'to-morrow and to-morrow '
subseq_ind = SubseqIndex(t, 8, 3)
occurrences, num_index_hits = querySubseqIndex(p,t,subseq_ind)
print(occurrences)
print(num_index_hits)  

t = open('1110.txt.utf-8').read()
p = 'English measure backward'
subseq_ind = SubseqIndex(t, 8, 3)
occurrences, num_index_hits = querySubseqIndex(p, t, subseq_ind)
print(occurrences)
print(num_index_hits)

# homework
t = readGenome('chr1.GRCh38.excerpt.fasta')

p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print('Q1: # of alignments are %d' % num_alignments)

p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print('Q2: # of character comparisons are %d' % num_character_comparisons)

p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
p_bm = BoyerMoore(p, alphabet='ACGT')
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print('Q3: # of alignments are %d' % num_alignments)

## p length 24, n = 2, k = 8, p = 'GGCGCGGTGGCTCACGCCTGTAAT'
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
all_matches,total_hits = approximate_matching_with_index(p, t, 2, 8)
print('Q4: # of occurrences with up to 2 mismatches are %d' % len(all_matches))

p = 'GGCGCGGTGGCTCACGCCTGTAAT'
all_matches,total_hits = approximate_matching_with_index(p, t, 2, 8)
print('Q5: # of total hits with up to 2 mismatches using pigeonhole approach are %d' % total_hits)


p = 'GGCGCGGTGGCTCACGCCTGTAAT'
subseq_ind = SubseqIndex(t, 8, 3)
occurrences, num_index_hits = querySubseqIndex(p,t,subseq_ind)
print('Q6: # of total hits with up to 2 mismatches using subseq is %d' % num_index_hits)  