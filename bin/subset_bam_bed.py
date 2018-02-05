#!/usr/bin/env python
"""
Subset a .bed file for only entries in chromosomes that correspond to chromosomes present in the given .bam files

workaround for issue here:
https://github.com/ding-lab/msisensor/issues/1

example usage:
./subset_bam_bed.py *.bam -b targets.bed -o targets_subset.bed
"""
import os
import sys
import argparse
import pysam


def get_bam_chroms(bam_file):
    """
    Gets all the chromosomes that a .bam file maps to
    """
    chroms = set()
    bam = pysam.AlignmentFile(bam_file, "rb")
    reads = bam.fetch()
    for read in reads:
        chroms.add(read.reference_name)
    return(chroms)

def subset_bed(bed_file, chroms, output_file):
    """
    Subset the .bed file for only lines that start with a chrom in the chroms list
    """
    with open(bed_file) as bed_in, open(output_file, "w") as bed_out:
        for line in bed_in:
            parts = line.strip().split('\t')
            if len(parts) > 0:
                chrom = parts[0]
                if chrom in chroms:
                    bed_out.write(line)

def main(bam_files, bed_file, output_file):
    """
    Main control function for the script
    """
    chroms = set()
    for bam_file in bam_files:
        bam_chroms = get_bam_chroms(bam_file)
        chroms.update(bam_chroms)
    subset_bed(bed_file = bed_file, output_file = output_file, chroms = bam_chroms)


def parse():
    """
    parse script args and pass them to `main`
    """
    # ~~~~ GET SCRIPT ARGS ~~~~~~ #
    parser = argparse.ArgumentParser(description = 'Script to subset msisensor microsatellites list file for only the chromosomes that are contained in the input bam files')

    # positional args
    parser.add_argument("bam_files", help = "Paths to input table files", nargs="+")

    # optional args
    parser.add_argument("-b", "--bed", nargs = 1,  required = True, dest = 'bed_file', metavar = 'bed_file', help="Path to the bed_file")
    parser.add_argument("-o", "--output", nargs = 1,  required = True, dest = 'output_file', metavar = 'output_file', help="name of output file")

    args = parser.parse_args()

    bam_files = args.bam_files
    bed_file = args.bed_file.pop()
    output_file = args.output_file.pop()
    main(bam_files, bed_file, output_file)


if __name__ == "__main__":
    parse()
