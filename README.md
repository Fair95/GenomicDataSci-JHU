# This is the Coursera course offered by JHU for Genomic Data Science


## Python in Genomic Data Science
### How to sort a dirtionary in Python
It is not possible to sort a dictionary, only to get a representation of a dictionary that is sorted. Dictionaries are inherently orderless, but other types, such as lists and tuples, are not. So you need an ordered data type to represent sorted values, which will be a list—probably a list of tuples.  

For instance,  

```python
import operator  
x = {1: 2, 3: 4, 4: 3, 2: 1, 0: 0}  
sorted_x = sorted(x.items(), key=operator.itemgetter(1))  
```
sorted_x will be a list of tuples sorted by the second element in each tuple. dict(sorted_x) == x.

And for those wishing to sort on keys instead of values:

```python
import operator
x = {1: 2, 3: 4, 4: 3, 2: 1, 0: 0}
sorted_x = sorted(x.items(), key=operator.itemgetter(0))
```

## Algorithms: Read Alignment Problem
### Naive matching
```python
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
```
Lets say y is the length of reference genome, x is the length of the reads:  
1. The maximum number of comparison would be **x(y-x+1)** times  
2. The minimum number of comparison would be **y-x+1** times  
3. This algorithm tends to have number of comparisons closer to minimum

### Boyer-Moore
Matching from left-to-right of the pattern and skip the maximum of either:  
1. Bad character rule: Upon mismatch, skip alignments until:  
  a) mismatch becomes a match, or  
  b) **p** moves past mismatched character
2. Good suffix rule: Let **t** = substring matched by inner loop; skip until:  
  a) there are no mismatches between **p** and **t** or  
  b) **p** moves past **t**  
3. naive matching (i.e. 1)

### Preprocessing
For a given pattern **p** and a given reference text **T**, pre-processing **p**, we can use the same pre-processed result for different **T**, samewise, we can use the same pre-processed result derived from **T** to test many different given **p**    
if an algorithm pre-processes **T**, then it is an ***Offline algorithm***  
Otherwise, it is an ***Online algorithm***, Boyer-Moore is an ***Online*** algorithm

### kmer indexing
k-mer: substring of length k  
Indexing DNA: Constructing an index which contain every possible k-mer substring in DNA **T** with its left-most offset associated.  
Quering index: Quering the any k-mer substring in pattern **p**, if we find one matching, we call it a hit which indicates there may be a matching around the region containing the substring.  
Verification: Chech the rest bases of the pattern matching or not.  
*Note kmer indexing is an offline algorithm (preprocessing **T**).*

### Distance
#### Hamming distance
For *X* & *Y* where |*X*| = |*Y*|,  
hamming distance = minimum # substitutions needed to turn one into the other.

#### Edit distance
For *X* & *Y*,  
edit distance = minimum # edits (substitutions, insertions, deletions) needed to turn one into the other

If len(*X*) != len(*Y*), the minimum of Edit distance is ||*X*| - |*Y*||, i.e. editDistance(*X*,*Y*) >= ||*X*| - |*Y*||

### Dynamic Programming for implementing Edit distance
Given a substring A and a character which in together comprise the string Aa, similarly, Bb is another string.  
edist(Aa,Bb) = min(edist(A,Bb)+1, edist(Aa,B)+1, edist(A,B)+delta(a,b))  
where delta = 1 if a != b, delta = 0 if a = b  

```python
def editDistance(x, y):
	D = []
	# initialize the matrix D with size (len(x)+1) * (len(y)+1)
	for i in range(len(x)+1):
		D.append([0]* (len(y)+1))

	# the distance between empty string and a string L is always len(L)
	for i in range(len(x)+1):
		D[i][0] = i
	for i in range(len(y)+1):
		D[0][i] = i

	for i in range(1, len(x)+1):
		for j in range(1, len(y)+1):
			distHor = D[i][j-1] + 1
			distVer = D[i-1][j] + 1
			if x[i-1] == y[j-1]:
				distDiag = D[i-1][j-1]
			else:
				distDiag = D[i-1][j-1] + 1

			D[i][j] = min(distHor, disVer, distDiag)
	return D[-1][-1]
```

#### Approximate matching in editDistance
For a [**p** X **T**] matrix:  
 we can find the minimal edit distance between **p** and a substring of **T** by  
Taking the minimum value in the bottom and then track back along with decreasing value.  
A vertical move in the traceback corresponds to a insertion in **p** with respect to **T**.  
A horizental move in tracback corresponds to a deletion in **p** with respect to **T**.

#### Global alignment --- Penalty matrix
Edit distance weighs all differences equally, where as global alignment can give different weights to different mutations
1. Transitions(A<->G or C<->T) are more likely to happen than transversions(others cases)  
2. Substitutions are more likely to happen than indels(deletions or insertions)

#### Local alignment --- scoring matrix
Global alignment is concerned with the overall similarity between two strings X and Y, whereas local alignment is concerned with similarities between substrings of X and Y
1. negative score for mismathes 
2. positive score for matches

## Algorithm: Assembly Problem
1. First law of assembly: If a suffix of read A is similar to a prefix of read B, then A and B might *overlap* in the genome.
2. Second law of assembly: More *coverage* leads to more and longer overlaps. *Coverage* = total length of all reads / length of the genome  
3. Third law of assembly: Reapeats makes Assembly Problem hard to solve since it is very hard to determine how many repeats are there.  
```python
def overlap(a, b, min_length = 3):
	start = 0

	while True:
		start = a.find(b[:min_length], start)
		if start == -1:
			return 0
		if b.startswith(a[start:]):
			return len(a)-start
		start += 1

from itertools import permutations

def naive_overlap_map(reads, k):
	olaps = {}
	for a,b in permutations(reads, 2):
		olen = overlap(a,b,min_length=k)
		if olen > 0:
			olaps[(a,b)] = olen
	return olaps
```

### Shortest Common Superstring
Give a set of string *S*, the *SCS* is the shortest string containing the strings in *S* as substrings
#### Directed Graph (Overlap Graph)
Each node represent a substring
Each arrow represent how many overlaps are between two strings(Directed: suffix of the previous overlap the prefix of the later)  

### Eulerian Walks
#### De Brujin Graphs --- multigraph
For every adjacent pair, adding an directed arrow.  
To use De Brujin Graph in Assembly problem, we:
1. First find all k_mers
2. Then find left k-1_mer and right k-1_mer in each k-mer
3. Lastly add directed edge to two k-1_mers

```python
def de_bruijn_ize(st, k):
	edges = []
	node = set()
	for i in range(len(st) - k + 1):
		edges.append((st[i:i+k-1]), st[i+1:i+k])
		nodes.add(st[i:i+k-1])
		nodes.add(st[i+1:i+k])
	return nodes, edges
```
Some properties of the graph:
1. Each k-mer has an edge
2. One node per distinct k-1_mer

To reconstruct the sequence of the assemble, just walk through the nodes along the directed edges exactly once. It is known as **Eulerian Walks**

## Problems in Real Assembly tools
1. Sequencing errors --- Deadend in graphs or disconnected graphs
2. Contain edges that can be inferred by other edges (transitional inferable) --- Occurs in overlap graph
3. Polyploidy --- Bubbles in the graphs
4. Repeats: Main problem, ususally assembles the unambiguous subsequences only.

## Possible Improvement
1. Increase the number of k in k-mers in order to contain the repeated patterns with unique patterns surrounded. (pair-end sequence)
