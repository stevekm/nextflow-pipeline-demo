#!/usr/bin/env perl


@fastqs = `find example-data/ -type f -name '*_R1_0*.fastq.gz' -or -name '*_R1.fastq.gz' -or -name '*_1.fastq.gz' | sort`;
@samples = ();

print "\@fastqs: ";
print @fastqs;

# print @fastqs;
while (my $fastq_r1 = shift(@fastqs)) {
    print "\n-----------------\n";
    # strip newlines from right end
    chomp($fastq_r1);

    print "\$fastq_r1: ";
    print "$fastq_r1\n";

    # print "\@fastqs: ";
    # print @fastqs;
    my $fastq_r2 = $fastq_r1;
    $fastq_r2 =~ s/(.*)_R1_0([0-9]+.fastq.gz)/${1}_R2_0${2}/;
    print "fastq_r2: $fastq_r2\n";
    $fastq_r2 =~ s/(.*)_R1.fastq.gz/${1}_R2.fastq.gz/;
    print "fastq_r2: $fastq_r2\n";
    $fastq_r2 =~ s/(.*)_1.fastq.gz/${1}_2.fastq.gz/;
    print "fastq_r2: $fastq_r2\n\n";


    my $sample = $fastq_r1;
    print "sample: $sample\n";
    # remove directory structure
    $sample =~ s/.*\///;
    print "sample: $sample\n";
    # bcl2fastq2 format (with S sample number)
    $sample =~ s/_S[0-9]{1,3}_L00[0-9]_R1.*//;
    print "sample: $sample\n";
    # bcl2fastq format with 2 barcodes
    $sample =~ s/_[ACTG]{6,}-[ACTG]{6,}_L00[0-9]_R1.*//;
    print "sample: $sample\n";
    # bcl2fastq format with 1 barcode
    $sample =~ s/_[ACTG]{4,}_L00[0-9]_R1.*//;
    print "sample: $sample\n";
    # no barcodes
    $sample =~ s/_L00[0-9]_R[12].*//;
    print "sample: $sample\n";
    # no barcodes or lane
    $sample =~ s/_R[12].fastq.gz//;
    print "sample: $sample\n";
    # no barcodes or lane
    $sample =~ s/_[12].fastq.gz//;
    print "sample: $sample\n";

    push @samples, $sample;

    print "\n-----------------\n";
}

print "samples: @samples\n";
#
# $bar = "foo bar";
# $bar =~ s/foo/bar/;
# print $bar;
