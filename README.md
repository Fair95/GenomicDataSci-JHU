## This is the Coursera course offered by JHU for Genomic Data Science

# How to sort a dirtionary in Python
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