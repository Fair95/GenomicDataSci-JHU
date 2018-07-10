## This script uses binary search package bisect to solve pattern find
## problem in genome.

import bisect

class Index(object):
	def __init__(self,t,k): # create the k-mer index object
							# t is the reference genome
							# k is the k_mer we defined
		self.k = k
		# initiate the k_mer index
		self.index = []

		for i in range(len(t) - k + 1):
			# the element in index is a tuple. e.g. ('GCT',3)
			# the first element in tuple is the k_mer, the second is
			# the left-most offset
			self.index.append((t[i:i+k], i))
		# sort the indext for downstreaming binary search
		self.index.sort()

	def query(self, p): # p is the pattern we need to find
		kmer = p[:self.k]
		i = bisect.bisect_left(self.index, (kmer,-1)) # -1 < 0 assures we always find the left most tuple
		hits = []
		# find all k_mers in indext after i that are the same as our pattern 
		while i < len(self.index):
			if self.index[i][0] != kmer:
				break
			hits.append(self.index[i][1])
			i += 1
		return hits
 
## Using subsequences instead of substring
class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
            	# print 'the reference is %s' %self.index[i][0]
            	# print 'the pattern is %s' % subseq
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def queryIndex(p,t,index):
	k = index.k
	offsets = []
	for i in index.query(p):
		# Verify if the rest of the bases match
		if p[k:] == t[i+k:i+len(p)]:
			offsets.append(i)
	return offsets

def HammingDistance(p, q):
    ham_dis = 0
    for i in range(len(p)):
        if (p[i] != q[i]):
            ham_dis += 1
    return ham_dis

def querySubseqIndex(p,t,index):
	k = index.k
	ival = index.ival
	offsets = set()
	total_hits = 0
	for i in range(ival):
		hits = index.query(p[i:])
		total_hits += len(hits)
		for j in index.query(p[i:]):
			if HammingDistance(p,t[j-i:j-i+len(p)])<=2:
				offsets.add(j-i)
	return list(offsets),total_hits