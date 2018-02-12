#!/usr/bin/env python

"""
Concatenates a list of files. The first line of the first file serves as a header for the output.

Usage:
$ find output/ -name "*.hg19_multianno.txt" | xargs ./concat_tables.py
"""
import argparse
import sys

def main(**kwargs):
    """
    Main control function for the script
    """
    file_list = kwargs.pop('file_list')
    output_file = kwargs.pop('output_file', None)

    if output_file:
        fout = open(output_file, "w")
    else:
        fout = sys.stdout

    # emit entire first file
    first_file = file_list.pop(0)
    with open(first_file) as fin:
        for line in fin:
            fout.write(line)

    # emit all lines except first from remaining files
    for input_file in file_list:
        with open(input_file) as fin:
            next(fin)
            for line in fin:
                fout.write(line)

    fout.close()


def parse():
    """
    parse script args and pass them to `main`
    """
    parser = argparse.ArgumentParser(description='This script will concatenate multiple table files with a common header')
    parser.add_argument("file_list", help="Paths to input table files", nargs="+")
    parser.add_argument("-o", default = False, dest = 'output_file', metavar = 'Table output file', help="Path to the output table file")
    args = parser.parse_args()

    main(**vars(args))

if __name__ == "__main__":
    parse()
