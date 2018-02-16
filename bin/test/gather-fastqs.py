#!/usr/bin/env python

"""
This script will read a list of samples from a samplesheet
and find all the fastq files in the supplied dir that match
and sort them into R1 and R2, if applicable
and output a new samplesheet with that information
"""
import argparse
import csv

def get_IEM_samples(input_sheet):
    """
    Get the sample IDs from the IEM format sheet
    https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf

    skip file contents down to this line in the samplesheet:

    [Data],,,,,,
    Sample_ID,Sample_Name,I7_Index_ID,index,Sample_Project,Description,GenomeFolder
    """
    samples = []

    # find the line that starts with "[Data]"
    data_start_line = 0
    with open(input_sheet) as f:
        for i, line in enumerate(f):
            parts = line.strip().split(',')
            print((i, parts))
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
            print(parts)


    return(samples)


def get_samples(input_sheet, sheet_format = "IEM"):
    """
    Get the samples from the samplesheet
    """
    samples = []
    if sheet_format == "IEM":
        samples = [sample for sample in get_IEM_samples(input_sheet)]
    return(samples)

def main(**kwargs):
    """
    Main control function for the script
    """

    print(kwargs)
    input_sheet = kwargs.pop('input_sheet')[0]
    output_file = kwargs.pop('output_file', "samples-fastq.csv")
    delim_input = kwargs.pop('delim_input', ",")
    delim_output = kwargs.pop('delim_output', ",")
    sheet_format = kwargs.pop('sheet_format', "IEM")

    sheet_samples = get_samples(input_sheet = input_sheet, sheet_format = sheet_format)
    output_samples = set(sheet_samples)

    # 

    print(sheet_samples)



def parse():
    """
    parse script args and pass them to `main`
    """
    parser = argparse.ArgumentParser(description='This script will find fastq files that match samples in the samplesheet')
    parser.add_argument("input_sheet", help="Paths to input samplesheet file", nargs=1)
    parser.add_argument("-o", default = "samples-fastq.csv", dest = 'output_file', metavar = 'output_file', help="Name of file to save new sheet to")
    parser.add_argument("-di", default = ",", dest = 'delim_input', metavar = 'delim_input', help="Input file delimiter")
    parser.add_argument("-do", default = ",", dest = 'delim_output', metavar = 'delim_output', help="Output file delimiter")
    parser.add_argument("-f", default = "IEM", dest = 'sheet_format', metavar = 'sheet_format', help="Format of input samplesheet")
    args = parser.parse_args()

    main(**vars(args))


if __name__ == '__main__':
    parse()
