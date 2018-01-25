params.bam_dir = "BAM-BWA"
params.pairs_sheet = "samples.pairs.csv"
params.regions_bed = "targets.bed"

Channel.fromPath( file(params.pairs_sheet) ) // read the .csv file into a channel
                        .splitCsv(header: true) // split .csv with header
                        .map { row -> // map each row in the .csv to files in the bam dir
                            sample_tumor_ID = row."${params.tumor_header}"
                            sample_normal_ID = row."${params.normal_header}"
                            return [
                            "sample_ID":"${sample_tumor_ID}_${sample_normal_ID}",
                            "sample_tumor_ID":"${sample_tumor_ID}",
                            "sample_tumor_bam":file("${params.bam_dir}/${sample_tumor_ID}.bam"),
                            "sample_tumor_bai":file("${params.bam_dir}/${sample_tumor_ID}.bam.bai"),
                            "sample_normal_ID":"${sample_normal_ID}",
                            "sample_normal_bam":file("${params.bam_dir}/${sample_normal_ID}.bam"),
                            "sample_normal_bai":file("${params.bam_dir}/${sample_normal_ID}.bam.bai")
                            ]
                        }
                        .into {sample_pairs; sample_pairs2}



process match_samples {
    executor "local"
    input:
    set val(sample_ID), val(sample_tumor_ID), file(sample_tumor_bam), file(sample_tumor_bai), val(sample_normal_ID), file(sample_normal_bam), file(sample_normal_bai) from sample_pairs

    exec:
    println "sample: ${sample_ID}, tumor: ${sample_tumor_ID}, sample_tumor_bam: ${sample_tumor_bam}"

}

process msisensor {
    clusterOptions '-pe threaded 1-32 -j y'
    publishDir "${params.output_dir}/${sample_ID}/MSI" , mode: 'copy', overwrite: true,
        saveAs: {filename ->
            if (filename == 'msisensor') "${sample_ID}.msisensor"
            else if (filename == 'msisensor_dis') "${sample_ID}.msisensor_dis"
            else if (filename == 'msisensor_germline') "${sample_ID}.msisensor_germline"
            else if (filename == 'msisensor_somatic') "${sample_ID}.msisensor_somatic"
            else null
        }

    input:
    file microsatellites from params.microsatellites
    set val(sample_ID), val(sample_tumor_ID), file(sample_tumor_bam), file(sample_tumor_bai), val(sample_normal_ID), file(sample_normal_bam), file(sample_normal_bai) from sample_pairs
    file regions_bed from file(params.regions_bed)

    output:
    file 'msisensor'
    file 'msisensor_dis'
    file 'msisensor_germline'
    file 'msisensor_somatic'

    script:
    """
    $params.msisensor_bin msi -d $microsatellites -n $bam_normal -t $bam_tumor -e $regions_bed -o msisensor -l 1 -q 1 -b \${NSLOTS:-1}
    """
}
