## Special Topics




*  Part 4: Special topics
    -  Biopython
        -  SeqIO, Maybe some of the NCBI support, converting files, 
            -  TODO: check the biopython tutorial/cookbook http://biopython.org/DIST/docs/tutorial/Tutorial.html
    -  CSV package maybe https://docs.python.org/3/library/csv.html



## Biopython

<img src="figures/biopython_logo_m.png" alt="BioPython" width="25%" align="right"/>


[Biopython](https://biopython.org/) is an open source collection of tools for manipulating biological data. The types of data it supports, and operations it provides for these data types are documented in the [Biopython Tutorial and Cookbook](https://biopython.org/DIST/docs/tutorial/Tutorial.html). 


### Getting started with BioPython

First start python (or ipython) if you haven't.
```bash
ll

cat helloworld.py

python3 helloworld.py
```

<div class="output">
shunter@tadpole:python$ ll
total 5.0K
drwxrwxr-x  2 shunter shunter 2.0K Oct  4 20:11 ./
drwxrwx--- 26 shunter shunter 2.0K Oct  4 20:14 ../
-rw-rw-r--  1 shunter shunter   23 Oct  4 20:09 helloworld.py<br>

shunter@tadpole:python$ cat helloworld.py 
print("Hello, World!")

shunter@tadpole:python$ python3 helloworld.py
Hello, World!

</div>


## Glob module