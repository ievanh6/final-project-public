#!/bin/bash

# trim paired end reads in parallel
# Info on the Trimmomatic tool available here:
# http://www.usadellab.org/cms/index.php?page=trimmomatic
# https://github.com/timflutre/trimmomatic

# Naupaka Zimmerman
# nzimmerman@usfca.edu
# November 6, 2017

# call script like so: "bash 06_trimmomatic.sh data/fastq/*_1.fastq"
# meant to run for paired reads based on name of forward
# read (ends in "_1.fastq"), which determines how variants of output are
# determined

# requires TruSeq3-SE.fa file with Illumina adapters to be in the directory
# this file is in the 'misc' folder for this project
# this set of parameters is a bit stricter than the defaults
# e.g. we're trimming any reads that drop below 15 in a window size of 4

for file in "$@"; do
	TrimmomaticSE -threads 8 $file \
		data/trimmed/$(basename -s _1.fastq $file).trim.fastq \
		ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 \
		SLIDINGWINDOW:4:15 MINLEN:36 &
done
