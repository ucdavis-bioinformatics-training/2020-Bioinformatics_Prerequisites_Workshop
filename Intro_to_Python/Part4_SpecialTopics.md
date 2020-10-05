## Special Topics




*  Part 4: Special topics
    -  Notes on documenting code
    -  Basic File Operations
    -  Biopython
        -  SeqIO, Maybe some of the NCBI support, converting files, 
            -  TODO: check the biopython tutorial/cookbook http://biopython.org/DIST/docs/tutorial/Tutorial.html
    -  CSV package maybe https://docs.python.org/3/library/csv.html



## Biopython

<img src="figures/biopython_logo_m.png" alt="BioPython" width="25%" align="right"/>


[Biopython](https://biopython.org/) is an open source collection of tools for manipulating biological data. The types of data it supports, and operations it provides for these data types are documented in the [Biopython Tutorial and Cookbook](https://biopython.org/DIST/docs/tutorial/Tutorial.html). 


### Getting started with BioPython

First start python (or ipython).

Lets start by exploring Biopython Sequence objects

```python
from Bio.Seq import Seq

m = Seq("AGTCCGTAG")

# Normal string indexing still work
print(m[0])

# Slicing still works, get the first two letters:
m[:2]

# And the last two:
m[-2:]

# Print the Positive and Negative number index
for index, letter in enumerate(m):
    print(f"{index} {index - len(m)} {letter}")

# In some cases sequences may need to be upper/lower case:
m.lower()

m.upper()

# Characters can be counted just like with strings:
100 * (m.count('C') + m.count('G')) / len(m)

# We can even join two Seq() objects, or a Seq() and a String:
n = m + "AAAAAAAAAAA"

print(n)

# However Seq objects also support a number of useful functions:
m.reverse_complement()


m.translate(table=11)


from Bio.Alphabet import IUPAC

n = Seq("ACGTCCGTCA", IUPAC.IUPACUnambiguousDNA)

```

<div class="output">

</div>


## Glob module


>>> from Bio.Seq import Seq
>>> from Bio.Alphabet import IUPAC
>>> my_seq = Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
...              IUPAC.protein)
>>> my_seq
Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF', IUPACProtein())
>>> print(my_seq)
MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
>>> my_seq.alphabet
IUPACProtein()
File:           /usr/lib/python3/dist-packages/Bio/Seq.py
Type:           type
Subclasses:     UnknownSeq
