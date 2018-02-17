#!/usr/bin/env python

"""
This script will read a list of samples from a samplesheet
and find all the fastq files in the supplied dir that match
and sort them into R1 and R2, if applicable
and output a new samplesheet with that information
"""
import argparse
import csv
import os
import sys
import find


def get_IEM_samples(input_sheet):
    """
    Get the sample IDs from the IEM format sheet
    https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf

    skip file contents down to this line in the samplesheet:

    [Data],,,,,,
    Sample_ID,Sample_Name,I7_Index_ID,index,Sample_Project,Description,GenomeFolder

    Parameters
    ----------
    input_sheet: str
        path to the IEM formatted samplesheet file

    Returns
    -------
    list
        a list of sample IDs from the sample sheet
    """
    samples = []

    # find the line that starts with "[Data]"
    data_start_line = 0
    with open(input_sheet) as f:
        for i, line in enumerate(f):
            parts = line.strip().split(',')
            if parts[0] == "[Data]":
                data_start_line = i
                break

    # start processing file from line after "[Data]"
    count = 0
    with open(input_sheet) as f:
        while count < data_start_line + 2:
            next(f)
            count += 1
        for line in f:
            parts = line.strip().split(',')
            sampleID = parts[0]
            samples.append(sampleID)
    return(samples)


def get_samples(input_sheet, sheet_format = "IEM"):
    """
    Get the samples from the samplesheet

    Parameters
    ----------
    input_sheet: str
        path to the samplesheet file
    sheet_format: str

    Returns
    -------
    list
        a list of sample IDs from the sample sheet
    """
    samples = []
    if sheet_format == "IEM":
        samples = [sample for sample in get_IEM_samples(input_sheet)]
    else:
        print("ERROR: unrecognized samplesheet format specified")
        raise

    return(samples)

def validate_samples_fastqs(samples_fastqs):
    """
    Checks for errors in the matching of the sample IDs and fastq files

    Parameters
    ----------
    samples_fastqs: list
        a list of dicts
    """
    # check to make sure no fastq file is listed for more than one sample
    # all samples in the list
    sample_IDs = [ sample_dict['sample'] for sample_dict in samples_fastqs ]
    for sample_ID in sample_IDs:
        # get the fastq files
        sample_all_fastqs = []
        for sample_dict in samples_fastqs:
            if sample_dict['sample'] == sample_ID:
                sample_all_fastqs = [ x for x in sample_dict['fastq-all'] ]
        print('{0}: {1}'.format(sample_ID, sample_all_fastqs))


def find_samples_fastqs(sample_IDs, fastq_dir, search_level = None):
    """
    Matches the supplied sample IDs to .fastq.gz files in the directory

    Parameters
    ----------
    sample_IDs: list
        a list of the ``str`` sample IDs to search for
    fastq_dir: str
        path to the directory to search for sample fastq.gz files
    search_level: int
        number of directories deep to search, or ``None``

    Returns
    -------
    list
        a list of dicts with information on the samples and their fastq files

    Notes
    -----
    Omits any supplied samples that did not have fastq files
    """
    samples_fastqs = []
    for sample_ID in sample_IDs:
        sample_pattern = '{0}*'.format(sample_ID)

        all_sample_fastqs = sorted(find.find(search_dir = fastq_dir, inclusion_patterns = [ sample_pattern, '*.fastq.gz' ], search_type = 'file', level_limit = search_level, match_mode = "all"))

        # dont return entries without fastq files
        if len(all_sample_fastqs) < 1:
            continue

        R1_fastqs = [ x for x in find.multi_filter(names = all_sample_fastqs, patterns = [ "*_R1_*.fastq.gz" ], match_mode = "all")]

        R2_fastqs = [ x for x in find.multi_filter(names = all_sample_fastqs, patterns = [ "*_R2_*.fastq.gz" ], match_mode = "all")]



        sample_dict = {'sample': sample_ID,
                        'fastq-all': all_sample_fastqs, # ','.join(all_sample_fastqs),
                        'fastq-R1': R1_fastqs, #','.join(R1_fastqs),
                        'fastq-R2': R2_fastqs #','.join(R2_fastqs)
                        }
        samples_fastqs.append(sample_dict)

    validate_samples_fastqs(samples_fastqs = samples_fastqs)

def main(**kwargs):
    """
    Main control function for the script
    """

    print(kwargs)
    input_sheet = kwargs.pop('input_sheet')[0]
    fastq_dir = kwargs.pop('fastq_dir')[0]
    output_file = kwargs.pop('output_file', "samples-fastq.csv")
    delim_input = kwargs.pop('delim_input', ",")
    delim_output = kwargs.pop('delim_output', ",")
    sheet_format = kwargs.pop('sheet_format', "IEM")
    search_level = kwargs.pop('search_level', None)
    # sample_IDs = kwargs.pop('sample_IDs', [])

    sheet_samples = get_samples(input_sheet = input_sheet, sheet_format = sheet_format)
    output_samples = set(sheet_samples)

    # make sure there were no duplicates by comparing lengths
    if len(sheet_samples) != len(output_samples):
        print("ERROR: The number of samples in the file does not match the number of unique samples. There may be duplicate samples in the input sheet")
        raise

    # get the fastq files
    find_samples_fastqs(sample_IDs = output_samples, fastq_dir = fastq_dir, search_level = search_level)





def parse():
    """
    parse script args and pass them to `main`
    """
    parser = argparse.ArgumentParser(description='This script will find fastq files that match samples in the samplesheet')
    parser.add_argument("input_sheet", help="Paths to input samplesheet file", nargs=1)
    parser.add_argument("fastq_dir", help="Paths to directory of fastq.gz files to search", nargs=1)
    parser.add_argument("-o", default = "samples-fastq.csv", dest = 'output_file', metavar = 'output_file', help="Name of file to save new sheet to")
    parser.add_argument("-di", default = ",", dest = 'delim_input', metavar = 'delim_input', help="Input file delimiter")
    parser.add_argument("-do", default = "\t", dest = 'delim_output', metavar = 'delim_output', help="Output file delimiter")
    parser.add_argument("--format", default = "IEM", dest = 'sheet_format', metavar = 'sheet_format', help="Format of input samplesheet")
    parser.add_argument("--search-level", default = None, dest = 'search_level', metavar = 'search_level', help="How many directories deep to search for matching files")
    # parser.add_argument("-s", default = [], dest = 'sample_IDs', metavar = 'sample_IDs', action='append', nargs = "*", help="Sample IDs to use instead of an input sheet file")

    args = parser.parse_args()

    main(**vars(args))


if __name__ == '__main__':
    parse()
