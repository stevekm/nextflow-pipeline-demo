#!/bin/bash

# USAGE: create_samplesheet.sh /path/to/fastq_dir
# creates a samplesheet to use for the pipeline

# script args
fastq_dir="${1:-none}"
[ "$fastq_dir" == 'none' ] && printf 'ERROR: fastq_dir not provided\n' && exit 1
[ ! -d "$fastq_dir" ] && printf "ERROR: fastq_dir is not a dir: ${fastq_dir}\n" && exit 1

[ -f "samples.fastq-raw.csv "] && rm -f "samples.fastq-raw.csv"
./gather-fastqs.pl "$fastq_dir"

sed 's/\,.*/,NA/g' "samples.fastq-raw.csv" | LC_ALL=C sort -u > "samples.tumor.normal.csv"
