def overlap(a, b, min_length = 3):
	""" Return length of longest suffix of 'a' matching
	 a prefix of 'b' that is at least 'min_length' characters long.
	 If no such overlap exists, return 0."""
	start = 0 # start all the way at the ledt

	while True:
		start = a.find(b[:min_length], start) # look for b's suffix in a
		if start == -1: # no more occurences to right
			return 0
		# found occurence; check for full suffix/prefix match
		if b.startswith(a[start:]):
			return len(a)-start
		start += 1 # move just past previous match

from itertools import permutations

def naive_overlap_map(reads, k):
	olaps = {}
	for a,b in permutations(reads, 2):
		olen = overlap(a,b,min_length=k)
		if olen > 0:
			olaps[(a,b)] = olen
	return olaps

def overlap_all_pairs(reads, k):
	kmer_dir = {}
	outgoing_node = set()
	for read in reads:
		for i in range(len(read) - k + 1):
			kmer_dir[read[i:i+k]] = set()
	for read in reads:
		for i in range(len(read) - k + 1):
			kmer_dir[read[i:i+k]].add(read)

	olaps = {}
	for a,b in permutations(reads,2):
		if b in kmer_dir[a[-k:]] and b != a:
			olen = overlap(a,b,min_length=k)
			if olen > 0:
				olaps[(a,b)] = olen
				outgoing_node.add(a)
	return olaps,outgoing_node

def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest

def pick_maximal_overlap(reads, k):
	reada, readb = None, None
	best_olen = 0
	for a, b in permutations(reads, 2):
		olen = overlap(a,b, min_length=k)
		if olen > best_olen:
			reada, readb = a, b
			best_olen = olen
	return reada, readb, best_olen

def pick_maximal_overlap_modified(reads, k):
	reada, readb = None, None
	best_olen = 0

	kmer_dir = {}
	outgoing_node = set()
	for read in reads:
		for i in range(len(read) - k + 1):
			kmer_dir[read[i:i+k]] = set()
	for read in reads:
		for i in range(len(read) - k + 1):
			kmer_dir[read[i:i+k]].add(read)

	for a, b in permutations(reads, 2):
		if b in kmer_dir[a[-k:]] and b != a:
			olen = overlap(a,b,min_length=k)
			if olen > best_olen:
				reada, readb = a, b
				best_olen = olen
	return reada, readb, best_olen
from copy import deepcopy
def pick_maximal_overlap_for_greedy_scs(reads, k, kmers):
	reada, readb = None, None
	best_olen = 0
	for read in reads:
		# only overlap the read_a with other reads that contain the suffix of the given read_a
		current_suffix = read[len(read)-k:]
		# retrive from k_mer dic for all reads containing suffix of read_a
		reads_with_kmer = deepcopy(kmers[current_suffix])
		# do not check overlap with itself
		reads_with_kmer.discard(read) 
		for read_with_kmer in reads_with_kmer:
			olen = overlap(read,read_with_kmer,k)
			if olen > best_olen:
				reada, readb = read, read_with_kmer
				best_olen = olen
	return reada, readb, best_olen
def greedy_scs_optimized(reads, k):
	kmers = {}
	# Constructing k_mer dic
	for read in reads:
		for i in range(len(read) - k + 1): # for each k-mer
			current_kmer = read[i:i+k]
			if kmers.has_key(current_kmer):
				kmers[current_kmer].add(read)
			else:
				kmers[current_kmer] = set([read])
	read_a, read_b, olen = pick_maximal_overlap_2(reads, k, kmers)
	while olen > 0:
		reads.remove(read_a)
		reads.remove(read_b)
		reads.append(read_a + read_b[olen:])
		# remove read_a and read_b from all k_mer dic
		for key, value in kmers.iteritems():
			for val in value.copy():
				if val == read_a or val == read_b:
					value.discard(val)
		read = read_a + read_b[olen:]
		# add assembled read to all k_mer dic it contains
		for i in range(len(read) - k + 1): # for each k-mer
			current_kmer = read[i:i+k]
			if kmers.has_key(current_kmer):
				kmers[current_kmer].add(read)
			else:
				kmers[current_kmer] = set([read])
		read_a, read_b, olen = pick_maximal_overlap_2(reads, k, kmers)
	return ''.join(reads)


def greedy_scs(reads, k):
	read_a, read_b, olen = pick_maximal_overlap(reads, k)
	while olen > 0:
		reads.remove(read_a)
		reads.remove(read_b)
		reads.append(read_a + read_b[olen:])
		read_a, read_b, olen = pick_maximal_overlap(reads,k)
	return ''.join(reads)

def de_bruijn_ize(st, k):
	edges = []
	node = set()
	for i in range(len(st) - k + 1):
		edges.append((st[i:i+k-1]), st[i+1:i+k])
		nodes.add(st[i:i+k-1])
		nodes.add(st[i+1:i+k])
	return nodes, edges

def scs_list(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    scs_trace = set()
    scs_set = set()
    for ssperm in permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup == None or len(sup) <= len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
            scs_trace.add(shortest_sup)
    for s in scs_trace:
    	if len(s) == len(shortest_sup):
    		scs_set.add(s)
    return scs_set  # return shortest