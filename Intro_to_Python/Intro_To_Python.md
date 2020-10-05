<img src="./figures/python-logo2x.png" alt="Python" width="40%" align="right"/>

<br><br><br>
------------
  
  
  
## Intro to Python



### Schedule
*  Part 1: Intro To Python and first program (1 - 1:30pm)
*  Part 2: Basics of Python and basic data types (1:30pm-3pm)
    *  Breakout session 1, Hands on exercises, Basic Data Types (1:45 - 2:00pm)
    *  Part 2 continued more advanced data types (2:00 - 2:15pm)
    *  Breakout session 2, Hands on exercises, Advanced Data Types (2:15 - 2:30pm)
    *  Part 2: Wrap up, objects and memory management 

3pm Break (15 minutes)

*  Part 3: Flow control (3:15 - 3:45pm) 
    *  Breakout session 3, Practice Flow Control (3:45 - 4pm)
*  Part 3: Continued, files, functions, documentation 
    *  Hands on Exercises 3
*  Part 4: Special topics (4 - 5pm)
    -  Biopython
    -  CSV package maybe https://docs.python.org/3/library/csv.html
    -  Plotting Seaborn / GGPlot





## Part 1: Intro to Python

### Why learn Python?


<img src="figures/programming_languages_recommended.png" alt="Python" width="70%" align="center"/>

https://businessoverbroadway.com/2019/01/13/programming-languages-most-used-and-recommended-by-data-scientists/


-  Python is extremely popular and widely used, especially for data science.
    +  Popular and getting more so in Bioinformatics, especially for building tools.
    +  For analysis, R (which you will learn later in the week) is arguably more useful currently due to the huge number of packages available from [Bioconductor](http://bioconductor.org/) and [CRAN](https://cran.r-project.org/). 
    +  The best option is to learn [Python, R, and bash](http://omgenomics.com/programming-languages/). A little of each will go a long way.
-  Freely available to [download](https://www.python.org/downloads/) for Windows, Linux, Mac OS X, etc.
-  Python is extremely versatile
    +  Used for a wide range of purposes from automating simple tasks to massive software projects with wide adoption [by many large companies.](https://realpython.com/world-class-companies-using-python/)
-  Installed on almost every Linux server.
-  Vast number of resources online: If you can Google for it you can learn how to do it.


------


### Goals for this section

-  Lower the barrier to entry by resolving basic "getting started" hurdles.
-  Provide tips and pointers on things to be aware of.
-  Provide a foundation in basic Python and some hands-on experience.
-  Give you some basic tools and recipes that you can build on in the future.
-  Get you hooked on programming and solving problems with Python.

-----

### Question
1)  Click "YES" if you have ever programmed in Python before.

2)  Click "YES" if you have ever programmed in any language before.

----- 

##  Background

###  What is a programming language and why do we need it?

Speaking to a computer in its native language is tedious and complicated. A programming language is a way for humans to describe a set of operations to a computer in a more abstract and understandable way. A helper program then translates our description of the operations into a set of instructions (machine code) for the computer to carry out.<br> 

Some day we may develop a programming language that allows us to communicate our instructions to the computer in our native language (Alexa, turn on the TV). Except for simple cases, this option doesn't exist yet, largely because human languages are complicated and instructions can be difficult to understand (even for other humans).


In order for the helper program to work properly, we need to use a concise language:
-  Well defined vocabulary for describing the basic set of supported operations.
-  Well defined set of Data Types that have a defined set of valid operations (add, subtract, etc).
-  Well defined syntax that leaves no ambiguity about when the computer should carry out each instruction.

Specifically in Python:

<img src="figures/python_interpreter.png" alt="PythonInterpreter" width="90%" align="center"/> <br><br>



-----

### A brief history of Python
-  Initially developed during the late 1980's by [Guido van Rossum](https://en.wikipedia.org/wiki/Guido_van_Rossum), BDFL until 2018.
-  First development version released in 1991. Version 1 released in 1994.
-  Python 2.0.0 released June, 2001
    -  Python 2.x end-of-life Jan 1, 2020.
    -  This version was so popular and widely used that many Bioinformatics programs were written using it. Some of these tools have been converted to support v3.x, others are in the process of being upgraded or have been abandoned and will stay on v2.x. The last Python 2.x release is still available for [download](https://www.python.org/downloads/release/python-2718/).
-  Python 3.x (December 2008) was a significant re-design and broke compatibility with some parts of v2.x.
-  The current version is 3.8.


------

### Interesting features of Python
-  High level: It hides a lot of the complicated details.
-  Interpreted: programs are compiled to byte code and run on a virtual machine instead of being compiled to native machine language
-  Garbage Collected: memory is allocated and freed for you automatically
-  Spaces matter in Python and are part of the language syntax. Be careful with copy/paste!
-  In Python, "Readability counts". 
    *  There is a style guide called [Python Enhancement Proposal 8 (PEP8)](https://www.python.org/dev/peps/pep-0008/) that documents and encourages consistent, readable code. If you plan to share your code with others, it is good practice to follow the style guide (or adopt the style used by the rest of the team). 
    *  These best practices are also known as writing "pythonic" or "idiomatic" python, [this guide](https://docs.python-guide.org/writing/style/) has more details. Try ```import this``` in your Python interpreter if you are a fan of programmer philosophy.

----

### Base Python and the extensive package ecosystem
-  Python has been extremely successful partly because it is modular and highly extensible. The core of Python is relatively small and compact, but this is supplemented with a large ["standard library"](https://docs.python.org/3/library/) that adds a large amount of additional functionality.
    +  Thousands of additional packages are available from the [PyPI](https://pypi.org/) repository.
    +  PythonPath variable
    +  Where do libraries live?
    +  Virtual Environments 
    +  Conflicts and package versions
        -  Virtual environments
        -  Conda



-  Key concepts of interacting with Python
    +  Interactive vs Script, python vs ipython





---------

### Writing Python Code

You can work directly on the interactive Python terminal, but eventually you might want a record of your work. One option is to edit your code in a text editor, then run it in the terminal. This is also a great way to test your code piece by piece as you work through developing or editing your program.


Lots of tools have been developed to help with writing computer code. Two popular options are [Visual Studio Code](https://code.visualstudio.com/), and [Atom Editor](https://atom.io/). Both work well and are available for free. Download and install either one of these if you want, or use a text editor you already have.

When working directly on the server, [Nano](https://www.nano-editor.org/), [Vim](https://www.vim.org/), and [Emacs](https://www.gnu.org/software/emacs/) are all popular editors.


-------


### Your first Python program

[Hello, World!](https://en.wikipedia.org/wiki/%22Hello,_World!%22_program) is traditionally the first program to write in any new programming language.


Connect to the server<br>
```bash
ssh <user>@tadpole.genomecenter.ucdavis.edu
```

Then follow this example:

```bash
export PATH=/share/biocore/shunter/opt/bin/:$PATH
alias ll='ls -alFh --color'

mkdir -p ~/python
cd ~/python

python3
```

<div class="output">
Python 3.8.2 (default, May  5 2020, 12:05:43) 
[GCC 5.4.0 20160609] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> 
</div>


Use the interactive Python interpreter to run a short program:
```python
>>> print("Hello World!")
Hello World!
```


<h3><font color="red">Once you have successfully run `python3` and printed "Hello World!", mark "Yes" in zoom. Post questions or problems to the Slack channel.</font></h3>

Now we need to exit the Python interpreter. This can be done by typing ```exit()``` or **CTRL+D**. 


Next, lets try writing a "Hello World" Python script.

```bash
nano helloworld.py
```

Enter the same print statement we used above, press Ctrl+O to write the changes, then Ctrl+X to exit the file.

Check that the file you created is present and that your changes were saved successfully, then run the script:

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

<h3><font color="red">Once you have successfully created your Python script and run it, mark "Yes" in zoom. Post questions or problems to the Slack channel.</font></h3>



----






------------- OLD STUFF for formatting examples --------------
### Intro

**Benchmarking Universal Single-Copy Orthologs**, [BUSCO](https://doi.org/10.1093/bioinformatics/btv351), is a popular software package for assessing genome/transcriptome assembly completeness using single copy orthologs. It was published in Oct 2015 and had 3486 citations as of July 2020 according to Google Scholar! [The authors](https://www.sib.swiss/evgeny-zdobnov-group) are also responsible for [OrthoDB](https://www.orthodb.org/), a large database of curated [orthologous genes](https://www.orthodb.org/orthodb_userguide.html#terminology). 

------


### How are BUSCOs made? <img src="figures/stork2.png" alt="busco_figure" width="30%" align="right"/>


The BUSCO sets are collections of nearly universally distributed (90%) single-copy orthologous genes found within species at a specific phylogenetic level. Originally these sets represented arthropods, vertebrates, metazoans, fungi, and eukaryotes, but additional genome sequences have made it possible to create BUSCO sets at a finer scale.

These sets are determined by analysis of species in the OrthoDB database. [The theory](https://academic.oup.com/gbe/article/doi/10.1093/gbe/evq083/573552) is that genes belonging to these sets are evolving under "single-copy control" where something about their necessity and dosage constraints maintains them at a single copy within the genome.

If a newly assembled genome or transcriptome is missing genes from the corresponding BUSCO set, something may have gone wrong with sequencing/assembly/annotation/etc, and other genes may also be missing.


1. Selection: Single-copy orthologs present in at least 90% of species in a specific group are selected from OrthoDB.

1. Model Building: Multiple sequence alignments of protein sequences from each BUSCO group are generated and used to build a hidden Markov model (HMM) for the group.

1. Pruning: Sequences are then searched against the library of HMM profiles to remove any that cannot reliably distinguish members of their group from other sequences.

1. Parameter Optimization: "expected-score" and "expected-length" classification cut-offs are calculated for each BUSCO group based on the distribution of HMM search scores and lengths for members of the group. These cut-offs well be used to classify new proteins as members of the group.

1. Consensus sequences and Block profiles (position-specific frequency matrices) for each BUSCO are then created.


-----

### How are genomes/transcriptomes assessed?

1. Consensus sequences for each BUSCO are searched against the genome sequence using tBLASTn. Regions containing potential BUSCOs are identified. Up to three candidate genomic regions can be identified for each BUSCO. 

1. Candidate regions are extracted from the genome and [AUGUSTUS/Augustus](http://augustus.gobics.de/) in combination with the BUSCO block profile is used for gene prediction. For transcriptomes,the protein prediction is used directly if available, otherwise the longest ORF within the transcript is used.

1. Each predicted gene is then matched against the BUSCO group's HMM profile, sequences meeting the minimum alignment cut-off are considered orthologous.

1. Orthologous sequences are then evaluated based on the expected-length cutoff. Sequences are classified as "Complete" if they meet the length cutoff, or "Fragmented" if too short. If multiple sequences meet the alignment and length cutoff they are classified as "Duplicated". Any BUSCO without Complete, Fragmented, or Duplicated sequence is "Missing".

1. Finally, "Complete" sequences are used to build a new gene prediction model for Augustus. A second round of Augustus gene prediction is then performed on all BUSCO-matching candidate regions that did not yield a "Complete" ortholog. Classification is then carried out a second time on the new set of predicted genes.

--------

BUSCO scores and contiguity as defined by N50 are not well correlated:

<img src="figures/Busco_vs_N50.png" alt="Busco_vs_N50" width="50%"/>

[Simão et al. 2015](https://academic.oup.com/bioinformatics/article/31/19/3210/211866)

---------

### Hands on, Activating/Installing BUSCO

First we need to set up BUSCO.

#### Create an interactive session:
```bash
srun -t 03:00:00 -c 20 -n 1 --mem 16000 --partition production --account genome_workshop --reservation genome_workshop --pty /bin/bash
aklog 
source ~/.bashrc  # only necessary if you have a ~/.bashrc

```

### Get access to BUSCO:

#### Option 1 Run BUSCO using modules 

```bash
mkdir -p /share/workshop/genome_assembly/$USER/busco
cd /share/workshop/genome_assembly/$USER/busco

module load busco/4.0.2

cp -r /share/biocore/shunter/2020-Genome_Assembly_Workshop/busco/augustus.config /share/workshop/genome_assembly/$USER/busco/
export AUGUSTUS_CONFIG_PATH=/share/workshop/genome_assembly/$USER/busco/augustus.config

cp /share/biocore/shunter/2020-Genome_Assembly_Workshop/busco/busco_config.ini /share/workshop/genome_assembly/$USER/busco/
export BUSCO_CONFIG_FILE=/share/workshop/genome_assembly/$USER/busco/busco_config.ini

cp /share/biocore/shunter/2020-Genome_Assembly_Workshop/busco/generate_plot.py /share/workshop/genome_assembly/$USER/busco/

busco --help

```

<font color="red">Once you have successfully run `busco --help` mark "Yes" in zoom. Post questions or problems to the Slack channel.</font>

-------

#### Option 2 Install BUSCO using Conda

<font color="red">This option is for patient people or people who need to install BUSCO on a system where no module is available.</font>

*If you go this route, you will not need to set environment variables or copy generate_plot.py as in Option 1.*

```bash
mkdir -p /share/workshop/genome_assembly/$USER/busco
cd /share/workshop/genome_assembly/$USER/busco
```

##### Download miniconda:

See: https://docs.conda.io/en/latest/miniconda.html for more details

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

##### Install it to your workshop folder:

```bash
sh Miniconda3-latest-Linux-x86_64.sh -b -p /share/workshop/genome_assembly/$USER/busco/miniconda
```

##### Activate this new Conda install:

```bash
eval "$(/share/workshop/genome_assembly/$USER/busco/miniconda/bin/conda shell.bash hook)"
```

##### Add some channels, update Conda:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda update --all
```

##### Create a new environment and install Busco:

Note that this step can take **a very long time** because Busco has a large number of dependencies. This set also sets the AUGUSTUS_CONFIG_PATH and BUSCO_CONFIG_FILE environment variables.

```bash
conda create -n busco_env
conda activate busco_env
conda install busco=4.0.6
```

----------

### Test BUSCO on a bacterial genome.

We need a genome to test on, lets start by assembling a small bacterial one.

The following code block symlinks in some raw Illumina reads and does some basic read clean up with [HTStream](https://github.com/s4hts/HTStream/issues).

```bash
# NOTE: Create an interactive session on the cluster if you closed the previous one.
cd /share/workshop/genome_assembly/$USER/busco
mkdir -p bacterial_test
cd bacterial_test

# Setup Raw Data
mkdir 00-RawData
ln -s /share/biocore/shunter/bacteria/*.gz ./00-RawData/

# Clean reads:
module load htstream/1.3.1
mkdir -p 01-HTS_Preproc

hts_Stats -L 01-HTS_Preproc/Bacteria.json -1 00-RawData/Bacteria_R1.fastq.gz -2 00-RawData/Bacteria_R2.fastq.gz | \
hts_SuperDeduper -A 01-HTS_Preproc/Bacteria.json | \
hts_SeqScreener -A 01-HTS_Preproc/Bacteria.json | \
hts_AdapterTrimmer -A 01-HTS_Preproc/Bacteria.json | \
hts_Stats -A 01-HTS_Preproc/Bacteria.json -F -f 01-HTS_Preproc/Bacteria

```

Next we assemble the cleaned reads with [Spades](http://cab.spbu.ru/software/spades/) and look at the assembly statistics.

```bash
module load spades/3.13.0
spades.py -t 20 -1 01-HTS_Preproc/Bacteria_R1.fastq.gz -2 01-HTS_Preproc/Bacteria_R2.fastq.gz -o 02-SpadesAssembly

module load assembly_stats/1.0.1

assembly_stats ./02-SpadesAssembly/contigs.fasta
```

```R
stats for ./02-SpadesAssembly/contigs.fasta
sum = 1113800, n = 60, ave = 18563.33, largest = 389847
N50 = 82313, n = 3
N60 = 75116, n = 5
N70 = 54132, n = 6
N80 = 37018, n = 9
N90 = 12161, n = 15
N100 = 128, n = 60
N_count = 0
Gaps = 0
```

<font color="red">Once you have successfully completed this step mark "Yes" in zoom. Post questions or problems to the Slack channel.</font>

-------

**Wow, an N50 of only 82Kb?**


#### Run BUSCO in genome assessment mode

We will use new features in BUSCO V4: better support for bacteria and archaea, auto-lineage selection, and automated download of reference datasets (all of which are very nice!). To speed things up we can ask Busco to only search the prokaryotic lineage using --auto-lineage-prok.

```bash 

busco -f -c 20 -m genome -i ./02-SpadesAssembly/contigs.fasta -o 03-Busco --auto-lineage-prok

```

```
        --------------------------------------------------
        |Results from dataset mycoplasmatales_odb10       |
        --------------------------------------------------
        |C:98.9%[S:98.3%,D:0.6%],F:1.1%,M:0.0%,n:174      |
        |172    Complete BUSCOs (C)                       |
        |171    Complete and single-copy BUSCOs (S)       |
        |1      Complete and duplicated BUSCOs (D)        |
        |2      Fragmented BUSCOs (F)                     |
        |0      Missing BUSCOs (M)                        |
        |174    Total BUSCO groups searched               |
        --------------------------------------------------
```

This isolate had previously been identified as *Mycoplasma ovipneumoniae* and Busco has identified it as part of the Mycoplasmatales family. The assembly looks like it captured almost all of the single copy genes. If we look into the Busco folders we can find some additional interesting information about the genome. Note that because this sample was a prokaryote Busco used, [Prodigal](https://github.com/hyattpd/Prodigal) to do gene prediction instead of Augustus.


Rather than use BUSCO's auto-lineage selection, we can also look through the BUSCO database and specify the lineage directly if we have a good identification for the sample. This will cause BUSCO to run more quickly:

```bash
busco --list-datasets

busco -f -c 20 -m genome -i ./02-SpadesAssembly/contigs.fasta -o 03-Busco_lineage --lineage_dataset mycoplasmatales_odb10

```

Finally we can generate the canonical BUSCO plot using a script that comes with the BUSCO package. However we need to install the ggplot2 package in R first.

Start R and run the following:

(answer yes to install the package to your personal library)

```R
install.packages("ggplot2")
q(save="no")
```

Next, copy the summary files and make the plot:

```bash
mkdir -p short_summaries
cp ./03-Busco/short_summary.* ./short_summaries/ 
cp ./03-Busco_lineage/short_summary.* ./short_summaries/
python3 /share/workshop/genome_assembly/$USER/busco/generate_plot.py -wd ./short_summaries/

```

<img src="figures/busco_figure.png" alt="busco_figure" width="80%"/>

Note that each of the summary files has been incorporated in the plot. This may be helpful in comparing different assemblies.

<font color="red">Once you have successfully completed this step mark "Yes" in zoom. Post questions or problems to the Slack channel.</font>

---------


## Team exercise: Test Busco on *Drosophila* HiFi assemblies.

-------

## <font color="red">First, exit your srun session if you haven't already. We will need all of the cluster resources for the next step.</font>

### Stop until `squeue --reservation genome_workshop` is empty

-------


Additional assemblies were built with:
1. Shasta version 0.5.1 using command:
    * shasta --input ELF_19kb.m64001_190914_015449.Q20.28X.fasta
1. Flye v2.7.1 using command:
    * python ./Flye/bin/flye -t 40 --pacbio-hifi ELF_19kb.m64001_190914_015449.Q20.28X.fasta --out-dir flyeasm

Setup: 

```bash 
cd /share/workshop/genome_assembly/$USER/busco/

mkdir -p drosophila_test
cd drosophila_test
```


#### Team 1: IPA primary contigs from Trio-binned Maternal assembly vs IPA Paternal assembly

HiFi reads for this assembly were binned based on Illumina kmers.

```
/share/workshop/genome_assembly/pacbio_2020_data_drosophila/hifi_long_read_mat_ipa_assembly/RUN/14-final/final.p_ctg.fasta
/share/workshop/genome_assembly/pacbio_2020_data_drosophila/hifi_long_read_pat_ipa_assembly/RUN/14-final/final.p_ctg.fasta
```

#### Team 2: IPA primary + accessory contigs vs Shasta contigs

```
/share/workshop/genome_assembly/pacbio_2020_data_drosophila/hifi_long_read_diploid_ipa_assembly/RUN/14-final/final.*.fasta
/share/biocore/shunter/drosophila/ShastaRun/Assembly.fasta
```


#### Team 3: IPA primary + accessory contigs vs Flye contigs

```
/share/workshop/genome_assembly/pacbio_2020_data_drosophila/hifi_long_read_diploid_ipa_assembly/RUN/14-final/final.*.fasta
/share/biocore/shunter/drosophila/flyeasm/assembly.fasta
```


#### Team 4: IPA primary contigs vs IPA accessory contigs after purge duplicates

```
/share/workshop/genome_assembly/pacbio_2020_data_drosophila/purge_dup_asm/final.purged.a_ctg.fasta
/share/workshop/genome_assembly/pacbio_2020_data_drosophila/purge_dup_asm/final.purged.p_ctg.fasta
```


#### Team 5: IPA contigs after purge_dups vs IPA primary contigs

```
/share/workshop/genome_assembly/pacbio_2020_data_drosophila/purge_dup_asm/final.purged.p_ctg.fasta
/share/workshop/genome_assembly/pacbio_2020_data_drosophila/hifi_long_read_diploid_ipa_assembly/RUN/14-final/final.p_ctg.fasta
```


#### Extra Credit:

```
/share/biocore/shunter/drosophila/ShastaRun/Assembly.fasta
/share/biocore/shunter/drosophila/flyeasm/assembly.fasta
```

-------------

### Team Instructions

### Part 1 (20 minutes in breakout rooms)

<font color="red">If problems or questions arise during this section, post them into the Slack channel and a TA will join your breakout room to help.</font>

1. Choose someone from the group to be the team lead and communicator. This person will be responsible for posting team results to the appropriate Slack thread, collating the results, and presenting them to the group. **This person should post a brief message to the Slack channel introducing their team and starting the team thread.**

1. Generate assembly statistics for the two genome assemblies you have been assigned. **Paste these into the Slack thread for your team. Make sure to label them clearly with the contig set!**

1. Discuss the assembly statistics as a group. Form some hypotheses about how the BUSCO scores will look. What are you expectations for a Drosophila genome? What do you know about how the assembly was done? **Record these hypotheses and your rational.**

1. As a team, figure out the proper commands for running BUSCO on the two assemblies assigned to your team. 
    * Use the "--lineage_dataset" option to speed up analysis. Figure out which BUSCO set is appropriate to use (hint: _Drosophila melanogaster_ is in order _Diptera_).
    * Configure BUSCO to use 40 cpus.
    * **Paste your solution into the Slack chat under your team's thread.**
    * Note: <font color="red">Don't start BUSCO yet (see next steps).</font>

1. Designate two people in from your team to run BUSCO. Assign one assembly to each.

1. Use srun to start an interactive session with **40 CPU** cores and 32 gigs of RAM. Start BUSCO and make sure that it is running correctly. (Alternatively build an sbatch script and submit the runs as a job). <font color="red">Remember to only submit two BUSCO jobs per team, there are not enough CPU resources for more than this.</font>

1. Announce in the Slack channel that your team has completed Part 1 and leave your breakout room.

-------

### Part 2 (10-20 minutes)

1. Once BUSCO has finished for both assemblies, paste the full path to the BUSCO results to your group's slack channel. Remember to label it clearly.

1. Aggregate the "short_summary.*" files from both of your contig sets. Make a BUSCO plot. Discuss and compare the results and evaluate your hypotheses.

1. Prepare a short presentation (~3 slides, 5 minutes) with your observations about the assembly statistics, initial hypotheses, BUSCO results and final conclusions.

1. Leave your breakout room and announce on the Slack channel that your group is ready to present.


### Extra Credit Part 3 (?? minutes)

1. Go through the Slack threads for each group. Collect that BUSCO output paths for each different set of contigs.

1. Collate all of the statistics and compare the results. Which assembler did better? 

<img src="figures/busco_figure_drosophila.png" alt="busco_figure" width="60%"/>

#### Assessment
**Flye contigs:**
* Very few missing or fragmented, some duplicates. Flye (in HiFi mode) is either discarding or collapsing many of the haplotigs.
    
**IPA_diploid_a:**
* These are the IPA accessory contigs before purging duplicates. Because IPA doesn't do a great job of binning haplotigs, the accessory set is missing many BUSCOs.

**IPA_diploid_p:**
* For the same reason that IPA accessory contigs are missing many BUSCOs, IPA primary contigs are excepted to have mostly complete BUSCOs but lots of duplicates. This is exactly what we see.
    
**IPA_diploid_a+p:**
* The combination of IPA_diploid_a + IPA_diploid_p (accessory and primary) contigs should result in fully complete, and mostly duplicated BUSCOs. This is exactly what we see.
    
**IPA_purged_a_ctg and IPA_purged_p_contig:**
* These sets of contigs have had purge duplicates run on them. This should remove duplicated sequce, resulting in sets of contigs that contain most BUSCOs in single copy.
* This is mostly true, however the accessory set is missing some, and also has a number of duplicates. Interestingly the total single copy BUSCOs (3206) + accessory duplicated BUSCOs sum to a larger number (3308) than the total number of BUSCOs in the set (3285).
    
**IPA_trio-mat and IPA_trio-pat:**
* For both of these contig sets, reads had been binned prior to assembly using the Trio binning approach. We would expect this binning to produce beautiful haploid assemblies. With some small exception, this is almost exactly what we see.
    
**Shasta:**
* Shasta is not designed for PacBio HiFi data and it shows. The Shasta assembly has the largest number of fragmented contigs, and is missinge a huge number of BUSCOs as well.