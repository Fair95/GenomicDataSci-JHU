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
	shortest_sup = None
	for ssperm in itertools.permutations(ss):
		sup = ssperm[0]
		for i in range(len(ss)-1):
			olen = overlap(ssperm[i], ssperm[i+1],min_length)
			sup += ssperm[i+1][olen:]
		if shortest_sup is None or len(sup) < len(shortest_sup):
			shortest_sup = sup
	return shortest_sup

def pick_maximal_overlap(reads, k):
	reada, readb = None, None
	best_olen = 0
	for a, b in itertools.permutations(reads, 2):
		olen = overlap(a,b, min_length=k)
		if olen > best_olen:
			reada, readb = a, b
			best_olen = olen
	return reada, readb, best_olen

def greedy_scs(reads, k):
	read_a, read_b, olen = pick_maximal_overlap(reads, k)
	while olen > 0:
		reads.remove(read_a)
		reads.remove(read_b)
		reads.append(read_a + read_b[olen:])
		read_a, read_b, olen = pick_maximal_overlap(reads,l)
	return ''.join(reads)

def de_bruijn_ize(st, k):
	edges = []
	node = set()
	for i in range(len(st) - k + 1):
		edges.append((st[i:i+k-1]), st[i+1:i+k])
		nodes.add(st[i:i+k-1])
		nodes.add(st[i+1:i+k])
	return nodes, edges