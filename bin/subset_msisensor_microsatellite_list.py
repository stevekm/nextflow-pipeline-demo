#!/usr/bin/env python
"""
Need to subset the msisensor microsatellites list for only contigs (chroms) that are present in the .bam input files
workaround for issue here:
https://github.com/ding-lab/msisensor/issues/1

"microsatellites.list" format:
chromosome      location        repeat_unit_length      repeat_unit_binary      repeat_times    left_flank_binary       right_flank_binary      repeat_unit_bases       left_flank_bases        right_flank_bases
chrM    65      1       2       6       439     953     G       CGTCT   TGTGC
chrM    202     2       11      3       448     783     GT      CTAAA   TAATT
chrM    208     4       195     3       955     638     TAAT    TGTGT   GCTTG

"""
import os
import sys
import argparse
import pysam
import csv

def index_bam(bam_file):
    """
    Runs samtools index on a .bam file
    """
    pysam.index(bam_file)

def bai_is_present(bam_file):
    """
    Checks to make sure the .bai for the given .bam file is present
    """
    bai_file = str(bam_file) + '.bai'
    bai_present = os.path.exists(bai_file)
    return(bai_present)

def validate_bam(bam_file):
    """
    Makes sure a .bam file is valid; put steps to validate a .bam file here
    """
    if not bai_is_present(bam_file):
        print('bai is not present for {0}, generating .bai index file'.format(bam_file))
        index_bam(bam_file)
        if not bai_is_present(bam_file):
            print('.bai file still not present, exiting')
            sys.exit(1)
        else:
            print('.bai generated')

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


def subset_microsatellites(microsatellites_file, output_file, bam_chroms):
    """
    Subset the microsatellites file for only the chromosomes that are in the .bam file
    """
    with open(microsatellites_file) as f_in, open(output_file, "w") as f_out:
        reader = csv.DictReader(f_in, delimiter = '\t')
        writer = csv.DictWriter(f_out, reader.fieldnames, delimiter = '\t')
        writer.writeheader()
        for row in reader:
            if row['chromosome'] in bam_chroms:
                writer.writerow(row)


def main(bam_files, microsatellites_file, output_file):
    """
    Main control function for the script
    """
    chroms = set()
    for bam_file in bam_files:
        validate_bam(bam_file)
        bam_chroms = get_bam_chroms(bam_file)
        chroms.update(bam_chroms)
    subset_microsatellites(microsatellites_file, output_file, bam_chroms)


def parse():
    """
    parse script args and pass them to `main`
    """
    # ~~~~ GET SCRIPT ARGS ~~~~~~ #
    parser = argparse.ArgumentParser(description='Script to subset msisensor microsatellites list file for only the chromosomes that are contained in the input bam files')

    # positional args
    parser.add_argument("bam_files", help="Paths to input table files", nargs="+")

    # optional args
    parser.add_argument("-m", "--microsatellites", nargs = 1,  required=True, dest = 'microsatellites_file', metavar = 'microsatellites_file', help="Path to the microsatellites_file")
    parser.add_argument("-o", "--output", nargs = 1,  required=True, dest = 'output_file', metavar = 'output_file', help="name of output file")

    args = parser.parse_args()

    bam_files = args.bam_files
    microsatellites_file = args.microsatellites_file.pop()
    output_file = args.output_file.pop()
    main(bam_files, microsatellites_file, output_file)


if __name__ == "__main__":
    parse()
