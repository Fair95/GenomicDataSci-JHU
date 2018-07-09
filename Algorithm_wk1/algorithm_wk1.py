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
			qual = f.readline.rstrip()
			# if we reach the end of the line,break and finish
			if len(seq) == 0:
				break
			sequences.append(seq)
			qualities.append(qual)
		return sequences, qualities

def Phred33ToQ(qual):
	return ord(qual) - 33

def QToPhred33(Q):
	return char(Q + 33)

def createHist(qualities):
	hist = [0] * 50 # lowest possible score is 2, highest possible score is 41, padding added just in case
	for qual in qualities:
		for phred in qual:
			q = Phred33ToQ(phred)
			hist[q] += 1
	import matplotlib.pyplot as plt
	plt.bar(range(len(hist)),hist)
	plt.show()

def findGCByPos(reads,length):
	gc = [0] * length
	totals = [0] * length
	for read in reads:
		len_read = len(read)
		for i in range(len_read):
			if read[i] == 'C' or read[i] == 'G':
				gc[i] += 1
			totals[i] += 1
	for i in range(len(gc)):
		if totals[i] > 0:
			gc[i] /= float(totals[i])
	import matplotlib.pyplot as plt
	plt.plot(range(len(gc)),gc)
	plt.show()
def countBases(seqs):
	import collections
	count = collections.Counter()
	for seq in seqs:
		count.update(seq)
	print (count)

def generateReads(genome, numReads, readLen):
	import random
	reads = []
	for _ in range(numReads):
		start = random.randint(0,len(genome)-readLen)
		reads.append(genome[start:start+readLen])
	return reads
