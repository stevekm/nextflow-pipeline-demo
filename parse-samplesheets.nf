params.samplesheet = "samples.pairs.controls.csv" // .csv samplesheet with tumor-normal sample ID pairs
params.bam_dir = "BAM-BWA" // directory containing .bam files for the samples


Channel.fromPath( file(params.samplesheet) ) // read the .csv file into a channel
                        .splitCsv(header: true) // split .csv with header
                        .map { row -> // map each row in the .csv to files in the bam dir
                            [
                            "sample_tumor_ID":"${row."#SAMPLE-T"}",
                            "sample_tumor_bam":file("${params.bam_dir}/${row."#SAMPLE-T"}.bam"),
                            "sample_tumor_bai":file("${params.bam_dir}/${row."#SAMPLE-T"}.bam.bai"),
                            "sample_normal_ID":"${row."#SAMPLE-N"}",
                            "sample_normal_bam":file("${params.bam_dir}/${row."#SAMPLE-N"}.bam"),
                            "sample_normal_bai":file("${params.bam_dir}/${row."#SAMPLE-N"}.bam.bai")
                            ]
                        }
                        .into{ sample_pairs; sample_pairs_subscr }

// testing out mapping .csv; print to console
sample_pairs_subscr.subscribe { row -> println "${row}"}


process match_samples {
    executor "local"
    input:
    set val(sample_tumor_ID), file(sample_tumor_bam), file(sample_tumor_bai), val(sample_normal_ID), file(sample_normal_bam), file(sample_normal_bai) from sample_pairs

    exec:
    println "tumor: ${sample_tumor_ID}, sample_tumor_bam: ${sample_tumor_bam}"

}
