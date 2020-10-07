A Simple Bioinformatics Workflow
=================================

**1\.** First, let's make a directory for our analysis:

	cd /share/workshop/prereq_workshop/$USER
	mkdir rnaseq
	cd rnaseq

Then, link in some raw data, take a look at the directory, and take a look at one of the files:

	ln -s /share/biocore/workshops/prereq_data/00-RawData
	ls -ltrh 00-RawData/
	zless 00-RawData/C61subset/C61subset_S67_L006_R1_001.fastq.gz

The raw data directory has two directories within it, each one corresponding to a sample. Within each of those is a subsetted set of paired-end reads for that sample.

---

**2\.** The first step in any analysis is preprocessing the data, such as removing adapters from the reads, trimming the reads based on quality, overlapping paired-end reads, etc. For this exercise we will just remove adapters from our reads. Make a directory for the trimmed output:

	mkdir 01-trimmed
	cd 01-trimmed

We will use a program from a suite called "htstream" to do the adapter trimming. Load the htstream module:

	module load htstream

Take a look at the options:

	hts_AdapterTrimmer -h

---

**3\.** Now, let's run it. We are going to use our subsetted raw input files, fastq output prefix, and a stats log file that gives us statistics about the state of the reads before and after the command. We are going to use the "--no-orphans" option to make things easier so that we only have to deal with a pair of files per sample.

	hts_AdapterTrimmer -F --read1-input ../00-RawData/C61subset/C61subset_S67_L006_R1_001.fastq.gz --read2-input ../00-RawData/C61subset/C61subset_S67_L006_R2_001.fastq.gz -fastq-output C61_trimmed --stats-file C61.stats.log

Take a look at the stats log file output:

	cat C61.stats.log

Now, go back to the command (by using the up arrow) and change it so that all the filenames are for the other sample. You can use <Ctrl>-Arrow keys to go left and right by entire words to make it easier to get to the parts you need to change. You need to figure this out and run it!

**\*\*FIGURE OUT THE COMMAND\*\***

Take a look at that stats file:

	cat I561.stats.log

---

**4\.** The next step is to do alignment, but instead of using a module, we are going to download and compile the software ourselves. So we are going to use an aligner called STAR, which is made for RNA-Seq data. Your challenge is to download and compile this software (in your software directory) and then add the path to STAR to your PATH variable in your .bash_profile file.

Here is the [STAR github page](https://github.com/alexdobin/STAR)

You will know if it worked if you can run STAR from any directory.

**\*\*YOU HAVE TO FIGURE OUT HOW TO INSTALL STAR\*\***

---

**5\.** Now we should be ready to run STAR. Go back to your rnaseq directory:

	cd /share/workshop/prereq_workshop/$USER/rnaseq

To run STAR we need an indexed genome to align against. So the first step is to get a genome, an annotation for that genome, and index it. Make a directory for the reference:

	mkdir ref
	cd ref

Link in reference and annotation files for Arabidopsis:

	ln -s /share/biocore/workshops/prereq_data/ref/* .

Make a directory for the STAR index files:

	mkdir star_index

Now take a look at the STAR options:

	STAR -h

There are a lot of them. In order to understand what is going on, you really need to read the manual for all of these software. So, go back to the github page for STAR and find the manual. Go to the section for generating indexes and look at the "Basic Options" section. We will use those options to generate the index:

	STAR --runMode genomeGenerate \
		--genomeDir star_index \
		--genomeFastaFiles Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
		--sjdbGTFfile Arabidopsis_thaliana.TAIR10.39.gtf \
		--sjdbOverhang 99 \
        --genomeSAindexNbases 12

This step takes about 8 minutes to run.

---

**7\.** Once we've generated an index, we are ready to align. Go back to your rnaseq directory:

	cd /share/workshop/prereq_workshop/$USER/rnaseq

And make a directory for the alignments:

	mkdir 02-alignment
	cd 02-alignment

Look at the manual section for "Running mapping jobs". We will be using those options, plus a few more. Run the command for both of the samples. Each one should take about a minute to run:

	STAR --genomeDir ../ref/star_index \
		--readFilesCommand zcat \
		--sjdbGTFtagExonParentGene gene_id \
		--sjdbGTFfile ../ref/Arabidopsis_thaliana.TAIR10.39.gtf \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode GeneCounts \
		--outFileNamePrefix C61_aligned \
		--readFilesIn ../01-trimmed/C61_trimmed_R1.fastq.gz ../01-trimmed/C61_trimmed_R2.fastq.gz

Now run the other one:

	STAR --genomeDir ../ref/star_index \
		--readFilesCommand zcat \
		--sjdbGTFtagExonParentGene gene_id \
		--sjdbGTFfile ../ref/Arabidopsis_thaliana.TAIR10.39.gtf \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode GeneCounts \
		--outFileNamePrefix I561_aligned \
		--readFilesIn ../01-trimmed/I561_trimmed_R1.fastq.gz ../01-trimmed/I561_trimmed_R2.fastq.gz

We will briefly discuss the output in class.
