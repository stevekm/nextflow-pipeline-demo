
params.samplesheet = "samples.pairs.controls.csv"

sample_pairs = Channel.fromPath( file(params.samplesheet) )
                        .splitCsv(header: true)
                        // .subscribe { row ->
                        //             println "${row."#SAMPLE-T"}\t${row."#SAMPLE-N"}"
                        //             // [${row."#SAMPLE-T"}, ${row."#SAMPLE-N"}]
                        //            }
println(sample_pairs)

// tuple = Channel.from( [1, 'alpha'], [2, 'beta'], [3, 'delta'] )
// println(tuple)

process match_samples {
    input:
    set val(sample_tumor), val(sample_normal) from sample_pairs

    script:
    """
    echo "tumor: ${sample_tumor}, normal: ${sample_normal}"
    """

}
