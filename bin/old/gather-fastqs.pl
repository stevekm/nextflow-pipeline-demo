#!/usr/bin/env perl

# script for parsing a directory of .fastq.gz files and generating a samplesheet
# borrowed from Igor; https://raw.githubusercontent.com/igordot/sns/master/gather-fastqs
# thanks Igor!

use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $HELP = <<HELP;

  Find FASTQ files in a given directory (must have "_R1" or "_1" in file name).
  Extract sample names and paired reads based on file names.
  Generate sample table file samples.fastq.csv in current directory.

  If run multiple times, it will add new samples to the sample table.

  usage: gather-fastqs dir

HELP

if (!$ARGV[0]) {
    die $HELP;
}

main();

# main subroutine
sub main {
    my $search_dir = $ARGV[0];

    # convert dir from relative to absolute
    $search_dir = `readlink -f $search_dir`;
    chomp($search_dir);

    # check that dir exists
    unless ( -d $search_dir ) {
        die "\n\n ERROR! $search_dir DOES NOT EXIST \n\n";
    }

    # find fastqs in given directory
    my $find_fastq_cmd_names = "-name '*_R1_0*.fastq.gz' -or -name '*_R1.fastq.gz' -or -name '*_1.fastq.gz'";
    my $find_fastq_cmd = "find -L $search_dir -maxdepth 3 -type f $find_fastq_cmd_names | LC_ALL=C sort";
    my @fastqs = `$find_fastq_cmd`;

    # counter single and paired reads
    my $reads_se = 0;
    my $reads_pe = 0;
    my @samples = ();

    # sample table file
    my $filename = "samples.fastq-raw.csv";
    open(my $fh, ">>", $filename);

    # process each fastq
    while (my $fastq_r1 = shift(@fastqs)) {
        chomp($fastq_r1);

        # check that R1 exists
        unless ( -e $fastq_r1 ) {
            die "\n\n ERROR! $fastq_r1 DOES NOT EXIST \n\n";
        }

        # generate R2 filename
        my $fastq_r2 = $fastq_r1;
        $fastq_r2 =~ s/(.*)_R1_0([0-9]+.fastq.gz)/${1}_R2_0${2}/;
        $fastq_r2 =~ s/(.*)_R1.fastq.gz/${1}_R2.fastq.gz/;
        $fastq_r2 =~ s/(.*)_1.fastq.gz/${1}_2.fastq.gz/;

        # blank if R2 does not exist
        unless ( -e $fastq_r2 ) {
            $fastq_r2 = "";
        }

        # blank if R2 is same as R1 (in case of not standard file name, for example)
        if ( $fastq_r1 eq $fastq_r2 ) {
            $fastq_r2 = "";
        }

        # count based on read type
        if ( length($fastq_r2) ) {
            $reads_pe++;
        }
        else {
            $reads_se++;
        }

        # extract sample name
        my $sample = $fastq_r1;
        # remove directory structure
        $sample =~ s/.*\///;
        # bcl2fastq2 format (with S sample number)
        $sample =~ s/_S[0-9]{1,3}_L00[0-9]_R1.*//;
        # bcl2fastq format with 2 barcodes
        $sample =~ s/_[ACTG]{6,}-[ACTG]{6,}_L00[0-9]_R1.*//;
        # bcl2fastq format with 1 barcode
        $sample =~ s/_[ACTG]{4,}_L00[0-9]_R1.*//;
        # no barcodes
        $sample =~ s/_L00[0-9]_R[12].*//;
        # no barcodes or lane
        $sample =~ s/_R[12].fastq.gz//;
        # no barcodes or lane
        $sample =~ s/_[12].fastq.gz//;

        push @samples, $sample;

        # show progress
        print STDERR " SAMPLE : $sample \n";
        print STDERR "  FASTQ R1 : $fastq_r1 \n";
        print STDERR "  FASTQ R2 : $fastq_r2 \n";

        # print sample table line
        my $output = "${sample},${fastq_r1},${fastq_r2}\n";
        print $fh "$output";

    }
    close($fh);

    # remove duplicate entries
    system("cat $filename | LC_ALL=C sort | uniq > ${filename}.tmp && mv -f ${filename}.tmp $filename");

    # get number of unique sample names
    my $num_files = @samples;
    @samples = uniq(@samples);
    my $num_samples = @samples;

    # print stats
    print STDERR "\n";
    print STDERR " NUMBER OF SAMPLES : $num_samples \n";
    print STDERR " NUMBER OF SINGLE FILES : $reads_se \n";
    print STDERR " NUMBER OF PAIRED FILES : $reads_pe \n";

}



# end
