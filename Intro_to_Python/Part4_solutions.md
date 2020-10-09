

### Exercises

### 1. Try to open the example_data.fastq.gz file without using gzip. Read the first line into a variable called "l", then print it. What happens?

```python
h = open("example_data.fastq.gz", 'rt')
l = next(h)
l
h.close()
```

### 2. How many base pairs long is the first read in the file?
```python
handle = gzip.open("example_data.fastq.gz", 'rt')
l = next(handle)
l = next(handle)

print(l)
print(len(l))  # 426bp!

handle.close()
```

### 3. The ord() function will return the ASCII value for a single character. Write some code to read the fourth line from the fastq file, then covert the values to a quality score and calculate the average quality value. 
    *  Hints: The last character is a line return (\n), use strip() to remove it. You can use list comprehension or a for loop to split up the characters and run ord() on each one.  

```python
handle = gzip.open("example_data.fastq.gz", 'rt')
l = next(handle)
l = next(handle)
l = next(handle)
l = next(handle)
handle.close()

print(l)

# Using list comprehension:
qvs = [ord(x)-33 for x in l.strip()]
sum(qvs)/len(qvs)

# Using for loops:
qvs = []
for ch in l.strip():
    qvs.append(ord(ch)-33)
sum(qvs)/len(qvs)

```
### 4. Write some code to figure out how many total lines are in example_data.fastq.gz. How many reads is this?

```python
# Write some code to figure out how many total lines are in example_data.fastq.gz. How many reads is this?
i = 0
for l in gzip.open("example_data.fastq.gz", 'rt'):
    i += 1
print(i)
```


----

### Exercises part2

### 1. Use SeqIO to filter the records in from example_data.fastq.gz. Create a new file called data_filtered.fastq. Only write reads that DO NOT start with "AGGG".

```python
from Bio import SeqIO
import gzip

with open("data_filtered.fastq", 'w') as outh:
    i = 0
    for record in SeqIO.parse(gzip.open("example_data.fastq.gz", 'rt'), 'fastq'):
        if record.seq[0:4] != 'AGGG':
            SeqIO.write(record, outh, 'fastq')
            i += 1
print(f"records: {i}")
```



### 2. Use SeqIO and the CSV package to record statistics for the reads in example_data.fastq.gz. Create a CSV file called read_stats.tsv. Create a DictWriter with column names: (ReadID, Length, GCcontent, AverageQuality). Collect each of these pieces of information for each read, and write it to the CSV file.  
```python
import csv
from Bio import SeqIO
import gzip

fieldnames = ['ReadID','Length','GCcontent', 'AverageQuality']
# Create a DictWriter, use tabs to delimit columns 
csvfile = open("read_stats.tsv", 'w')
writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
writer.writeheader()

for record in SeqIO.parse(gzip.open("example_data.fastq.gz", 'rt'), 'fastq'):
    data={}
    data['ReadID'] = record.name
    data['Length'] = len(record)
    data['GCcontent'] = 100 * (record.seq.count('C') + record.seq.count('G')) / len(record)
    data['AverageQuality'] = sum(record.letter_annotations['phred_quality'])/len(record) 
    writer.writerow(data)

csvfile.close()

```

### 3. Use SeqIO and a Dictionary, count the frequency of the first 15bp and last 15bp of each read in example_data.fastq.gz. What are the most common sequences?


```python
from Bio import SeqIO
import gzip
import pprint

# Create dictionary objects to store data in
first15 = {}
last15 = {}
first_and_last15 = {}

# Iterate over records, extract first and last 15bp, count them
for record in SeqIO.parse(gzip.open("example_data.fastq.gz", 'rt'), 'fastq'):
    f15 = str(record.seq[0:15])
    l15 = str(record.seq[-15:])

    if f15 not in first15:
        first15[f15] = 0
    first15[f15] += 1

    if l15 not in last15:
        last15[l15] = 0
    last15[l15] += 1

    if (f15,l15) not in first_and_last15:
        first_and_last15[(f15,l15)] = 0
    first_and_last15[(f15,l15)] += 1

# Sort and print with pprint (pretty printer). Note that pprint actually sorts dictionaries
pprint.pprint(sorted(first15.items(), key=lambda x: x[1], reverse=True))

pprint.pprint(sorted(last15.items(), key=lambda x: x[1], reverse=True))

pprint.pprint(sorted(first_and_last15.items(), key=lambda x: x[1], reverse=True))


```