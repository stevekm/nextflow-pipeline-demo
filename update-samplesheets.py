#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script will update the samplesheets ``samples.analysis.json`` and ``samples.analysis.tsv`` output by the ``generate-samplesheets.py`` script.
Use this script to update those files with tumor-normal pairs sample metadata, from the file ``samples.tumor.normal.csv`` (output by Excel)


"""
NA_value = "NA"
delim = ','
tumor_normal_sheet = 'samples.tumor.normal.csv'
samples_analysis_json = 'samples.analysis.json'
samples_analysis_tsv = 'samples.analysis.tsv'

tumor_normal_samples = []

with open(tumor_normal_sheet) as f:
    for line in f:
        parts = line.strip().split(delim)
        print(parts)
        if len(parts) >= 2:
            tumor_ID = parts[0]
            normal_id = parts[1]
            line_dict = {}

            print()
