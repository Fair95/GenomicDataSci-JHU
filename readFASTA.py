import sys
import string
import getopt
import operator
def usage():
	print """
	This should gives the helper messge!
	"""

o, a = getopt.getopt(sys.argv[1:],'l:f:h') # o = list of optional arguments
										 # a = list of required arguments
										 #l: means l is an optinal argument requiring value
										 #h means h is an optional argument not requiring value
opts = {}
frame = -1
seqlen = 0

# add the optinal arguments user provided to the dictionary (easier to retrive)
for k,v in o:
	opts[k] = v
# h means helper message
if '-h' in opts.keys():
	usage(); sys.exit()
# we need at least one required argument (which is the filename)
if len(a) < 1:
	usage(); sys.exit("input fasta file is missing")
else:
	filename = a[0]
# if user provide optional arguments l (length of the sequence), we set the seqlen correspondingly
if '-l' in opts.keys():
	if int(opts['-l'])<0:
		print "Length of sequence should be positive!"; sys.exit(0)
	seqlen = opts['-l']
## if user provide optional arguments f(frame), we set the frame accordingly
if '-f' in opts.keys():
	if int(opts['-f'])<=0 or int(opts['-f'])>3:
		print "Reading frame should be either 1 or 2 or 3"; sys.exit(0)
	frame = opts['-f']

try:
	f = open(filename,'r')
except IOError:
	print "File %s does not exist!!" % filename
seqs = {}
for line in f:
	# discard the newline at the end (if any)
	line = line.rstrip()
	# distinguish header from sequence
	if line[0] == '>': # or line.startswith('>')
		words = line.split()
		# discard '>' sign
		name = words[0][1:]
		seqs[name] =''
	else: # sequence, not header
		seqs[name] = seqs[name]+line
f.close()


print  "number of records are", len(seqs.keys())
len_rank = []
for key in seqs:
	print "sequence with identifier", key, "has the a length of ", len(seqs[key])
	len_rank.append(len(seqs[key]))
print "the length of sequences are listed from lowest to highest:",sorted(len_rank)
def find_stop_codon(dna,search_pos):
	stop_codon_pos = -1
	stop_codon = ['TAA','TAG','TGA']
	for i in range(search_pos,len(dna),3):
		codon = dna[i:i+3].upper()
		if codon in stop_codon:
			stop_codon_pos = i
			break
	return stop_codon_pos

def find_start_codon(dna,search_pos):
	start_codon_pos = -1
	start_codon = ['ATG']
	for i in range(search_pos,len(dna),3):
		codon = dna[i:i+3].upper()
		if codon in start_codon:
			start_codon_pos = i
			break
	return start_codon_pos

def find_codon(dna,search_pos):
	start_pos = find_start_codon(dna,search_pos)
	stop_pos = find_stop_codon(dna,start_pos)
	if start_pos == -1:
		return [-1,-1]
	elif stop_pos == -1:
		return [-1,-1]
	else:
		return [start_pos,stop_pos]

## The followling section aims to find all ORF (from start codon to end codon) in a given dna sequence
if int(frame)!= -1:
	ORF_rank = []
	for key in seqs:
		search_pos = int(frame)-1
		print "============================"
		while find_codon(seqs[key],search_pos)[0]!= -1:
			start_pos = find_codon(seqs[key],search_pos)[0]
			stop_pos = find_codon(seqs[key],search_pos)[1]
			print "sequence with identifier", key,"has ORF", seqs[key][start_pos:stop_pos+3], "with length of", stop_pos+3-start_pos, "starting at position", start_pos+1
			ORF_rank.append(stop_pos+3-start_pos)
			search_pos = stop_pos+3
	print sorted(ORF_rank,reverse = 1)[0]
if int(seqlen)!=0:
	dna_dir = {}
	for key in seqs:		
		for i in range(len(seqs[key])-int(seqlen)):
			repeat = seqs[key][i:i+int(seqlen)]
			if repeat not in dna_dir:
				dna_dir[repeat] = 1
			else:
				dna_dir[repeat] += 1
	print sorted(dna_dir.items(), key=operator.itemgetter(1),reverse=1)


