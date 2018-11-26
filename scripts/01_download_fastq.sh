#!/bin/bash

# the pipe and tail -n +2 is a handy way to exclude the first line
for SRA_number in $(cut -f 6 data/metadata/SraRunTable.txt | tail -n +2)
do
    fastq-dump.2.9.2 -v $SRA_number -O data/raw_data
done
