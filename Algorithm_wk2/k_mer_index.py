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

def queryIndex(p,t,index):
	k = index.k
	offsets = []
	for i in indext.query(p):
		# Verify if the rest of the bases match
		if p[k:] == t[i+k:len(p)]:
			offsets.append(i)
	return offsets
