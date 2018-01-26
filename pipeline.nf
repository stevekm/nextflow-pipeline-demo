params.input_dir = "sns-dir"
params.regions_bed = "targets.bed"
params.probes_bed = "probes.bed"

params.bam_bwa_dir = "${params.input_dir}/BAM-BWA"
params.pairs_sheet = "${params.input_dir}/samples.pairs.csv"
params.variants_sheet = "${params.input_dir}/summary.VCF-GATK-HC-annot.csv"
params.variants_gatk_hc_dir = "${params.input_dir}/VCF-GATK-HC"

//
// read sample pairs from samplesheet
//
Channel.fromPath( file(params.pairs_sheet) ) // read the .csv file into a channel
                        .splitCsv(header: true) // split .csv with header
                        .map { row -> // map each row in the .csv to files in the bam dir
                            sample_tumor_ID = row."${params.tumor_header}"
                            sample_normal_ID = row."${params.normal_header}"
                            return [
                            "${sample_tumor_ID}_${sample_normal_ID}", // sample comparison ID
                            "${sample_tumor_ID}",
                            file("${params.bam_bwa_dir}/${sample_tumor_ID}.bam"),
                            file("${params.bam_bwa_dir}/${sample_tumor_ID}.bam.bai"),
                            "${sample_normal_ID}",
                            file("${params.bam_bwa_dir}/${sample_normal_ID}.bam"),
                            file("${params.bam_bwa_dir}/${sample_normal_ID}.bam.bai")
                            ]
                        }
                        .into {sample_pairs_demo;
                            sample_pairs_msi}

//
// read samples from variants summary sheet
//
Channel.fromPath( file(params.variants_sheet) )
                    .splitCsv(header: true)
                    .map{row ->
                        sample_ID = row."${params.variant_sample_header}"
                        return [
                        sample_ID,
                        file("${params.variants_gatk_hc_dir}/${sample_ID}.vcf")
                        ]
                    }
                    .into { sample_variants }




//
// pipeline steps
//
process match_samples {
    executor "local"
    input:
    set val(sample_ID), val(sample_tumor_ID), file(sample_tumor_bam), file(sample_tumor_bai), val(sample_normal_ID), file(sample_normal_bam), file(sample_normal_bai) from sample_pairs_demo

    exec:
    println "sample: ${sample_ID}, tumor: ${sample_tumor_ID}, sample_tumor_bam: ${sample_tumor_bam}"

}

process msisensor {
    clusterOptions '-pe threaded 1-8 -j y'
    publishDir "${params.output_dir}/MSI" , mode: 'move', overwrite: true, //"${params.output_dir}/${sample_ID}/MSI"
        saveAs: {filename ->
            if (filename == 'msisensor') "${sample_ID}.msisensor"
            else if (filename == 'msisensor_dis') "${sample_ID}.msisensor_dis"
            else if (filename == 'msisensor_germline') "${sample_ID}.msisensor_germline"
            else if (filename == 'msisensor_somatic') "${sample_ID}.msisensor_somatic"
            else null
        }

    input:
    file microsatellites from file(params.microsatellites)
    set val(sample_ID), val(sample_tumor_ID), file(sample_tumor_bam), file(sample_tumor_bai), val(sample_normal_ID), file(sample_normal_bam), file(sample_normal_bai) from sample_pairs_msi
    file regions_bed name "regions.bed" from file(params.regions_bed)

    output:
    file 'msisensor'
    file 'msisensor_dis'
    file 'msisensor_germline'
    file 'msisensor_somatic'

    script:
    """
    $params.msisensor_bin msi -d $microsatellites -n $sample_normal_bam -t $sample_tumor_bam -e $regions_bed -o msisensor -l 1 -q 1 -b \${NSLOTS:-1}
    """
}

process deconstructSigs_signatures {
    publishDir "${params.output_dir}/genomic_profiles" , mode: 'move', overwrite: true,
        saveAs: {filename ->
            if (filename == 'sample_signatures.Rds') "${sample_ID}_signatures.Rds"
            else if (filename == 'sample_signatures.pdf') "${sample_ID}_signatures.pdf"
            else if (filename == 'sample_signatures_pie.pdf') "${sample_ID}_signatures_pie.pdf"
            else null
        }
    input:
    set val(sample_ID), file(sample_vcf) from sample_variants

    output:
    file 'sample_signatures.Rds'
    file 'sample_signatures.pdf'
    file 'sample_signatures_pie.pdf'

    script:
    """
    $params.deconstructSigs_make_signatures_script "$sample_vcf"
    """
}
