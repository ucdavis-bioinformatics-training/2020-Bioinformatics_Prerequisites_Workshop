## Special Topics


## Outline 

*  [Part 4: Special topics](#top)
    -  [Basic File Operations](#files)
    -  [Biopython](#Biopython)
    -  [CSV package](#CSV)
 


## <a name="top">Basic File Operations</a>

Reading and writing files in Python is easy with built in support for common operations.

In order to work with a file, we send a request to the operating system to open a **handle** to the file. This **handle** allows us to send and/or receive information from the file depending on what **mode** we open the file in. 

The steps are:
1) Open a file handle.
2) Read or write data.
3) Close file handle.

#### Writing to a file:

```python
# 1) Open a file handle to write text to:
handle = open("test_file.txt", "w")

# 2) Lets practice looping over a range of integers and writing them to our file:
for i in range(0,20):
    handle.write(f"Integer {i}\n")

# 3) Close the file handle
handle.close()

# If you are running ipython try "cat test_file.txt" to check the contents of the file.
```

#### Reading from a file:

Now lets reverse the process and read the information back in from the file. This time we will try an alternative approach that automatically closes the file for us when we are done:

```python

# Option 1, iterates over the file line by line:
with open("test_file.txt", "r") as handle:
    for line in handle:
        print(line.strip())

# Option 2, reads the entire file in and prints it
with open("test_file.txt", "r") as handle:
    print(handle.read())

# Bonus: What if we want to look at the file manually?
handle = open("test_file.txt", "r")
line = next(handle)
print(line)
handle.close()

```

-------

#### What about compressed files?

```python
# First we need a compressed file to work with.
# Use urllib to download the file from GitHub

import urllib
url = "https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2020-Bioinformatics_Prerequisites_Workshop/master/Intro_to_Python/example_data.fastq.gz"

urllib.request.urlretrieve(url, 'example_data.fastq.gz')


# Use more functions in os to check on the download:
import os
os.path.isfile("example_data.fastq.gz")

os.path.getsize("example_data.fastq.gz")

# If you are running ipython, try "less -S example_data.fastq.gz" to get a look at the file.
```

The file we just download is a [FastQ](https://en.wikipedia.org/wiki/FASTQ_format) formatted file, containing data from an Illumina MiSeq sequencer. If you look closely, you will see that each record is made up of four lines. You will learn more about this file type later in the workshop.

1) The first line starts with "@" and contains the Read ID.
2) The second line contains the actual sequence data (ACGT or N).
3) The third line is a "+" sign, and in some older files also contains the Read ID.
4) The fourth line contains an [ASCII encoding of the Phred + 33 quality score](https://en.wikipedia.org/wiki/FASTQ_format#Quality). 

```
@M02034:478:000000000-CYBKG:1:1101:26272:4914_P3:Z:seq1 1:N:0:61
AGACAGGCAGTGGGAGCAAGATGGATTCACAGGCCCAGGTTCTTATATTGCTGCTGCTATGGGTATCTGGTACCTGTGGGGACATTGTGATGTCACAGTCTCCATCCTCCCTGGCTGTGTCAGCAGGAGAGATGGT
+
:FDGGGCGGDGEDEEEGGGGGGGGGGGGGGGGFGGEFGGEFFGGFGFDFFF@FGFGGCABFF?GGCFGGGCFCGGFGGGG@EGGGGFGDGFGGFGGGGGF<AFGGGGGGGGGFGF@@EGGGCFFGCFCGGGGGGGF
@M02034:478:000000000-CYBKG:1:1101:4726:7224_P3:Z:seq1 2:N:0:61
AGAGATGGAGACAGACACACTCCTGTTATGGGTACTGCTGCTCTGGGTTCCAGGTTCCACTGGTGACATTGTGCTGACACAGTCTCCTGCTTCCTTAGCTGTATCTCTGGGGCAGAGGGCCACCATCTCATACAGG
+
GGGGGGGGGGGGGGGGGGGGGGGGGGGGFGFFGGGGFGGGGGCFGGGGGGGGGGGGGGGGGGGGGDEGGGFFG<EFGGGF@GGGGGGGD<FFF=FDFFGFFGGGGFEFFGGGGGGGGGGGIIIII"IICIIIIIII
@M02034:478:000000000-CYBKG:1:1101:27026:8436_P3:Z:seq1 2:N:0:61
GAAGCATCCTCTCTTCCAGCTCTCAGAGATGGAGACAGACACACTCCTGTTATGGGTACTGCTGCTCTGGGTTCCAGGTTCCACTGGTGACATTGTGCTGACACAGTCTCCTGCTTCCTTAGCTGTATCTCTGGGG
+
GGGGGGGGGGGGGGGGGGGG<FEFGGGC8BEFCF8FFGEFGGGEF?FACG99@DFBFGGG<FAEFCFGGGFFFGGEC<C=FGFGFGGGGGGGFFAEGGGGGGGGGG?EFGFGGGGFGGGGGGGGGFCEE,FFGDGG
```


If we try to open this file like we did with the previous file we will get incomprehensible garbage. 

```python
# EXERCISE: Try to open the file, read one line into a variable called "l" and then print l. What happens?


```

In order to read the file properly, we need to uncompress it. The Python gzip module allows us to do this.


```python
import gzip

# Note that 'rt' opens the file for text reading
handle = gzip.open("example_data.fastq.gz", 'rt')

l = next(handle)

print(l)

```


Exercises

1) How many base pairs long is the first read in the file?

2) The ord() function will return the ASCII value for a single character. Can you read the fourth line from the fastq file and calculate the average quality value? Hints: The last character is a line return, use strip() to remove it. You can use list comprehension or a for loop to split up the characters and run ord() on each one.  

3) Write some code to figure out how many total lines are in example_data.fastq.gz. How many reads is this?

4) 

```python



```

## <a name="top">BioPython</a>

<img src="figures/biopython_logo_m.png" alt="BioPython" width="25%" align="right"/>


[Biopython](https://biopython.org/) is an open source collection of tools for manipulating biological data. The types of data it supports and operations it provides for these data types are documented in the [Biopython Tutorial and Cookbook](https://biopython.org/DIST/docs/tutorial/Tutorial.html) and are quite extensive. 


Explore Biopython Sequence objects

```python
from Bio.Seq import Seq

# Create a Seq object:
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


```


[SeqIO](https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec51) provides an interface to read/write a large number of different Bioinformatics file formats. 

Using SeqIO should look very familiar after learning about file handles. Note that we need to provide a file type so that SeqIO knows how to properly read our file.

```python 
# Import SeqIO module
from Bio import SeqIO

# Create  a gzip file handle 
handle = gzip.open("example_data.fastq.gz", 'rt')

# Tell SeqIO to parse the data in the file in fastq format
parser = SeqIO.parse(handle, 'fastq')

# Extract the first record
record = next(parser)

# Explore the record object:
record.id

record.seq

len(record)

# Calculate average quality score
sum(record.letter_annotations['phred_quality'])/len(record) 


```

Use SeqIO to get the lenths of all sequences in the file and calculate statistics.

```python
from Bio import SeqIO
import gzip
import statistics

with gzip.open("example_data.fastq.gz", 'rt') as h:
    l = [len(r) for r in SeqIO.parse(h, 'fastq')]


max(l)

min(l)

sum(l)/len(l)

# How many reads are < 427bp?
sum(m < 429 for m in l)

sum(m > 429 for m in l)

```

Another common use of SeqIO is to convert from one file type to another:

```python
from Bio import SeqIO
import gzip

with gzip.open("example_data.fastq.gz", 'rt') as inf, gzip.open("example_data.fasta.gz", 'wt') as outf:
    SeqIO.convert(inf, "fastq", outf, "fasta")

# If you are running ipython, use less to check whether the file has been converted.
```


---------

## <a name="top">CSV package</a>

It is commonly necessary to read or write tabular data when doing bioinformatics work. Python has a handy module to support this called [CSV](https://docs.python.org/3/library/csv.html).

Lets experiment with writing and reading tabular files using the DictWriter and DictReader classes. 

Writing data out.
```python
import csv

# First lets write an example file
with open('my_data.csv', 'w') as csvfile:
    # Create a list with the column headers
    fieldnames = ['SampleID','Age','Treatment', 'Weight']
    # Create a DictWriter, use tabs to delimit columns 
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
    # Write out the header
    writer.writeheader()
    
    writer.writerow({"SampleID":"mouse1", "Age":4, "Treatment":"Control", "Weight":3.2})
    writer.writerow({"SampleID":"mouse2", "Age":5, "Treatment":"Control", "Weight":3.6 })
    writer.writerow({"SampleID":"mouse3", "Age":4, "Treatment":"Control", "Weight":3.8 })
    writer.writerow({"SampleID":"mouse4", "Age":4, "Treatment":"ad libitum", "Weight":3.6 })
    writer.writerow({"SampleID":"mouse5", "Age":4, "Treatment":"ad libitum", "Weight":3.7 })
    writer.writerow({"SampleID":"mouse6", "Age":4, "Treatment":"ad libitum", "Weight":3.5 })

# If you are running ipython, try checking the contents of the file with "cat my_data.csv"

```

Reading data in from a file
```python
import csv

with open('my_data.csv', 'r') as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    for row in reader:
        print(f"SampleID:{row['SampleID']}, Age:{row['Age']}, Treatment:{row['Treatment']}, Weight:{row['Weight']}")


```


Now lets try a little 

```python
# Create an empty dictionary
chrcount = {}

# read in each line of the fasta file, count character occurrences
with gzip.open('example_data.fasta.gz', 'rt') as h:
    for line in h:
        for ch in line.strip():
            if ch not in chrcount:
                chrcount[ch] = 0
            chrcount[ch] += 1

# Print out the character (dictionary key) and associated value (the count):
for k, v in chrcount.items():
    print(f"character:{k}, count:{v}")

```

What if we want to print out values from most common to least? First we need a way to sort, we can use Python's "sorted" function. The sorted() function is a little bit tricky because the "key" argument has to be a function. This is very powerful because it allows flexibility, but can also be confusing.

```python
# Check how sorted works:
sorted([4,2,5,6,2])

# 
sorted([4,2,5,6,2], reverse=True)

# Sorting your dictionary by value requires some tricks:
sorted(chrcount.items(), key=lambda x: x[1], reverse=True)

# How does this work?
def getv(x):
    return(x[1])

sorted(chrcount.items(), key=getv, reverse=True)


```


Exercises

1) Use SeqIO and the CSV package to record statistics for the reads in example_data.fastq.gz. Create a CSV with column names, (ReadID, Length, GCcontent, ) 

2) Use SeqIO and a Dictionary, count the frequency of the first 15bp and last 15bp of each read in example_data.fastq.gz.   