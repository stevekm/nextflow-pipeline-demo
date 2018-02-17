#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script will search the provided directories for fastq.gz files matching naming criteria
and generate sample sheets to use for analysis, in .tsv and .json format.
Multiple directories can be supplied, and only the unique files will be used.

Script overview:
- search for all .fastq.gz files

Usage
-----
Example usage:

    ./gather-fastqs2.py example-data example-data example-data


Notes
------
Old samplesheets with the same name will be overwritten

Only file names as described below are supported:

    https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm
    Naming Convention

    FASTQ files are named with the sample name and the sample number, which is a numeric assignment based on the order that the sample is listed in the sample sheet. For example: Data\Intensities\BaseCalls\SampleName_S1_L001_R1_001.fastq.gz
    	▶ 	SampleName—The sample name provided in the sample sheet. If a sample name is not provided, the file name includes the sample ID, which is a required field in the sample sheet and must be unique.
    	▶ 	S1—The sample number based on the order that samples are listed in the sample sheet starting with 1. In this example, S1 indicates that this sample is the first sample listed in the sample sheet.

    Note

    Reads that cannot be assigned to any sample are written to a FASTQ file for sample number 0, and excluded from downstream analysis.
    	▶ 	L001—The lane number.
    	▶ 	R1—The read. In this example, R1 means Read 1. For a paired-end run, there is at least one file with R2 in the file name for Read 2. When generated, index reads are I1 or I2.
    	▶ 	001—The last segment is always 001.
"""

import os
import sys
import find
import re
import csv
import json


# ~~~~~ FUNCTIONS ~~~~~ #
def collapse(dicts, collapse_key):
    """
    Collapses a list of dicts based on a given key

    Parameters
    ----------
    dicts: list
        a list of dicts that are assumed to have common keys
    collapse_key: str
        dict key to collapse on

    Returns
    -------
    list
        a list of dicts

    Examples
    --------
    Example usage

        x = [{'ID': 'foo', 'File1': 'foo-1-1.txt', 'File2': 'foo-1-2.txt'},
        {'ID': 'foo', 'File1': 'foo-2-1.txt', 'File2': 'foo-2-2.txt'},
        {'ID': 'foo', 'File1': 'foo-3-1.txt', 'File2': 'foo-3-2.txt'},
        {'ID': 'bar', 'File1': 'bar-1-1.txt', 'File2': 'bar-1-2.txt'},
        {'ID': 'bar', 'File1': 'bar-2-1.txt', 'File2': 'bar-2-2.txt'},
        {'ID': 'bar', 'File1': 'bar-3-1.txt', 'File2': 'bar-3-2.txt'}]
        collapse(x, 'ID')

        # {'File2': ['foo-1-2.txt', 'foo-2-2.txt', 'foo-3-2.txt'], 'File1': ['foo-1-1.txt', 'foo-2-1.txt', 'foo-3-1.txt'], 'ID': 'foo'}
        # {'File2': ['bar-1-2.txt', 'bar-2-2.txt', 'bar-3-2.txt'], 'File1': ['bar-1-1.txt', 'bar-2-1.txt', 'bar-3-1.txt'], 'ID': 'bar'}
    """
    ids = list(set([i[collapse_key] for i in dicts])) # unique values from all dicts for the given key
    keys = set([k for k in i.keys() for i in dicts]) # all keys in all dicts
    collapsed_list = [] # list to hold output dicts
    for i in ids:
        y = {collapse_key: i} # retain value of collapsed field
        for k in keys :
            if k != collapse_key: # add the rest of the keys
                if k not in y.keys():
                    y.update({k: []}) # initialize them as lists for appending
                for z in dicts: # iterate over the original list of dicts again
                    if z[collapse_key] == i: # only the entries that match the selected collapse value
                        y[k].append(z[k]) # append value to the list
        collapsed_list.append(y) # add the collapsed dict to the list
    return(collapsed_list)


def main():
    """
    Main control function for the script
    """
    # get script args
    search_dirs = sys.argv[1:]
    fastqs_R1 = []
    samples = []

    if len(search_dirs) < 1:
        print("ERROR: no directories provided")
        sys.quit()

    # find the R1 fastq files
    for search_dir in search_dirs:
        for fastq_R1 in sorted(find.find(search_dir = search_dir,
                                        inclusion_patterns = [ '*_R1_0*.fastq.gz' ],
                                        search_type = 'file' )):
            fastqs_R1.append(fastq_R1)
    fastqs_R1 = list(sorted(set((fastqs_R1))))

    for R1_name in fastqs_R1:
        # generate R2 filename
        R2_name = os.path.join(os.path.dirname(R1_name), re.sub(r'(.*)_R1_0([0-9]+\.fastq\.gz)', r'\1_R2_0\2', os.path.basename(R1_name)))
        if not os.path.exists(R2_name): R2_name = None

        # extract sample name
        sample_name = re.sub(r'_S[0-9]{1,3}_L00[0-9]_R1.*', '', os.path.basename(R1_name))
        sample_dict = {
            'Sample': sample_name,
            'R1': R1_name,
            'R2': R2_name
        }
        samples.append(sample_dict)

    # save long version of the table; one line per R1 R2 pair
    with open('samples.fastq-long.tsv', 'w') as f:
        writer = csv.DictWriter(f, delimiter= '\t', fieldnames=['Sample', 'R1', 'R2'])
        writer.writeheader()
        for item in samples:
            writer.writerow(item)

    # save a JSON
    with open('samples.fastq-long.json', 'w') as f:
        json.dump(samples, f, sort_keys = True, indent = 4)


    # reduce to condensed version; one entry per sample with all R1 and R2
    samples_collapsed = collapse(dicts = samples, collapse_key = 'Sample')

    # add some extra metadata
    for sample_dict in samples_collapsed:
        sample_dict['Tumor'] = sample_dict['Sample']
        sample_dict['Normal'] = None

    # save a JSON
    with open('samples.analysis.json', 'w') as f:
        json.dump(samples_collapsed, f, sort_keys = True, indent = 4)

    # prepare dicts for .tsv printing
    samples_to_print = []
    for sample_dict in samples_collapsed:
        # convert fastq lists to comma separated
        new_R1 = ','.join([str(x) for x in sample_dict['R1']])
        new_R2 = ','.join([str(x) for x in sample_dict['R2']])
        sample_dict['R1'] = new_R1
        sample_dict['R2'] = new_R2
        samples_to_print.append(sample_dict)

    # write to file
    with open('samples.analysis.tsv', 'w') as f:
        writer = csv.DictWriter(f, delimiter= '\t', fieldnames=['Sample', 'Tumor', 'Normal', 'R1', 'R2'])
        writer.writeheader()
        for item in samples_to_print:
            writer.writerow(item)


if __name__ == '__main__':
    main()
