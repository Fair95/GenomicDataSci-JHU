# This is the Coursera course offered by JHU for Genomic Data Science

<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

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

## Algorithms in Genomic Data Science
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

If len(*X*) != len(*Y*), the minimum of Edit distance is ||*X*| - |*Y*||, i.e. $$editDistance(*X*,*Y*) \ge ||*X*| - |*Y*||$$

