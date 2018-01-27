params.input_dir = "sns-dir"
params.regions_bed = "targets.bed"
params.probes_bed = "probes.bed"

params.bam_bwa_dir = "${params.input_dir}/BAM-BWA"
params.bam_gatk_ra_rc_dir = "${params.input_dir}/BAM-GATK-RA-RC"

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
                            file("${params.bam_gatk_ra_rc_dir}/${sample_tumor_ID}.dd.ra.rc.bam"),
                            "${sample_normal_ID}",
                            file("${params.bam_gatk_ra_rc_dir}/${sample_normal_ID}.dd.ra.rc.bam"),
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
                    .into { sample_variants_deconstructSigs }

Channel.fromPath( file(params.variants_sheet) )
                    .splitCsv(header: true)
                    .map{row ->
                        sample_ID = row."${params.variant_sample_header}"
                        // println [sample_ID, row].join('\t')
                        println "-------------"
                        println "--- start csv row mapping ----"
                        println "row: " + row
                        println 'sample_ID: ' + sample_ID
                        println "${params.bam_gatk_ra_rc_dir}/${sample_ID}.dd.ra.rc.bam"
                        println file("${params.bam_gatk_ra_rc_dir}/${sample_ID}.dd.ra.rc.bam")
                        println "--- end csv row mapping ----"
                        println "-------------"

                        return [
                        sample_ID,
                        file("${params.bam_gatk_ra_rc_dir}/${sample_ID}.dd.ra.rc.bam")                        ]
                    }
                    .into {
                        sample_bam_delly2;
                        sample_bam_delly2_deletions;
                        sample_bam_delly2_duplications;
                        sample_bam_delly2_inversions;
                        sample_bam_delly2_translocations;
                        sample_bam_delly2_insertions;
                        sample_bam_gatk_coverage_custom;
                        sample_bam_demo
                    }

//
// pipeline steps
//
process match_samples {
    tag { sample_ID }
    executor "local"
    input:
    set val(sample_ID), val(sample_tumor_ID), file(sample_tumor_bam), val(sample_normal_ID), file(sample_normal_bam) from sample_pairs_demo

    exec:
    println "sample: ${sample_ID}, tumor: ${sample_tumor_ID}, sample_tumor_bam: ${sample_tumor_bam}"

}

process check_samples_mapping {
    tag { sample_ID }
    executor "local"
    input:
    set val(sample_ID), file(bam_gatk_ra_rc) from sample_bam_demo

    exec:
    println "sample: ${sample_ID}, bam: ${bam_gatk_ra_rc}"
}


process msisensor {
    tag { sample_ID }
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
    set val(sample_ID), val(sample_tumor_ID), file(sample_tumor_bam), val(sample_normal_ID), file(sample_normal_bam) from sample_pairs_msi
    file regions_bed name "regions.bed" from file(params.regions_bed)

    output:
    file 'msisensor'
    file 'msisensor_dis'
    file 'msisensor_germline'
    file 'msisensor_somatic'

    script:
    """
    $params.samtools_bin index ${sample_normal_bam}
    $params.samtools_bin index ${sample_tumor_bam}
    $params.msisensor_bin msi -d $microsatellites -n $sample_normal_bam -t $sample_tumor_bam -e $regions_bed -o msisensor -l 1 -q 1 -b \${NSLOTS:-1}
    rm -f ${sample_normal_bam}.bai
    rm -f ${sample_tumor_bam}.bai
    """
}

process deconstructSigs_signatures {
    tag { sample_ID }
    errorStrategy 'ignore' // script fails if there are no variants
    publishDir "${params.output_dir}/Genomic_Profiles" , mode: 'move', overwrite: true,
        saveAs: {filename ->
            if (filename == 'sample_signatures.Rds') "${sample_ID}_signatures.Rds"
            else if (filename == 'sample_signatures.pdf') "${sample_ID}_signatures.pdf"
            else if (filename == 'sample_signatures_pie.pdf') "${sample_ID}_signatures_pie.pdf"
            else null
        }
    input:
    set val(sample_ID), file(sample_vcf) from sample_variants_deconstructSigs

    output:
    file 'sample_signatures.Rds'
    file 'sample_signatures.pdf'
    file 'sample_signatures_pie.pdf'

    script:
    """
    $params.deconstructSigs_make_signatures_script "$sample_vcf"
    """
}

// process delly2 {
//     tag { sample_ID }
//     publishDir "${params.output_dir}/SNV-Delly2", mode: 'move', overwrite: true,
//         saveAs: {filename ->
//             if (filename == 'deletions.vcf') "${sample_ID}_deletions.vcf"
//             else if (filename == 'duplications.vcf') "${sample_ID}_duplications.vcf"
//             else if (filename == 'inversions.vcf') "${sample_ID}_inversions.vcf"
//             else if (filename == 'translocations.vcf') "${sample_ID}_translocations.vcf"
//             else if (filename == 'insertions.vcf') "${sample_ID}_insertions.vcf"
//             else null
//         }
//
//     input:
//     set val(sample_ID), file(bam_gatk_ra_rc), file(bai_gatk_ra_rc) from sample_bam_delly2
//
//     output:
//     file "deletions.vcf"
//     file "duplications.vcf"
//     file "inversions.vcf"
//     file "translocations.vcf"
//     file "insertions.vcf"
//
//     script:
//     """
//     $params.delly2_bin call -t DEL -g ${params.hg19_fa} -o deletions.bcf "${bam_gatk_ra_rc}"
//     $params.delly2_bcftools_bin view deletions.bcf > deletions.vcf
//
//     $params.delly2_bin call -t DUP -g ${params.hg19_fa} -o duplications.bcf "${bam_gatk_ra_rc}"
//     $params.delly2_bcftools_bin view duplications.bcf > duplications.vcf
//
//     $params.delly2_bin call -t INV -g ${params.hg19_fa} -o inversions.bcf "${bam_gatk_ra_rc}"
//     $params.delly2_bcftools_bin view inversions.bcf > inversions.vcf
//
//     $params.delly2_bin call -t BND -g ${params.hg19_fa} -o translocations.bcf "${bam_gatk_ra_rc}"
//     $params.delly2_bcftools_bin view translocations.bcf > translocations.vcf
//
//     $params.delly2_bin call -t INS -g ${params.hg19_fa} -o insertions.bcf "${bam_gatk_ra_rc}"
//     $params.delly2_bcftools_bin view insertions.bcf > insertions.vcf
//     """
// }


process delly2_deletions {
    tag { sample_ID }
    publishDir "${params.output_dir}/SNV-Delly2-deletions", mode: 'move', overwrite: true,
        saveAs: { filename -> "${sample_ID}_deletions.vcf" }

    input:
    set val(sample_ID), file(bam_gatk_ra_rc) from sample_bam_delly2_deletions

    output:
    file "deletions.vcf"

    script:
    """
    $params.samtools_bin index ${bam_gatk_ra_rc}
    $params.delly2_bin call -t DEL -g ${params.hg19_fa} -o deletions.bcf "${bam_gatk_ra_rc}"
    $params.delly2_bcftools_bin view deletions.bcf > deletions.vcf
    rm -f ${bam_gatk_ra_rc}.bai
    """
}

process delly2_duplications {
    tag { sample_ID }
    publishDir "${params.output_dir}/SNV-Delly2-duplications", mode: 'move', overwrite: true,
        saveAs: { filename -> "${sample_ID}_duplications.vcf" }

    input:
    set val(sample_ID), file(bam_gatk_ra_rc) from sample_bam_delly2_duplications

    output:
    file "duplications.vcf"

    script:
    """
    $params.samtools_bin index ${bam_gatk_ra_rc}
    $params.delly2_bin call -t DUP -g ${params.hg19_fa} -o duplications.bcf "${bam_gatk_ra_rc}"
    $params.delly2_bcftools_bin view duplications.bcf > duplications.vcf
    rm -f ${bam_gatk_ra_rc}.bai
    """
}

process delly2_inversions {
    tag { sample_ID }
    publishDir "${params.output_dir}/SNV-Delly2-inversions", mode: 'move', overwrite: true,
        saveAs: { filename -> "${sample_ID}_inversions.vcf" }

    input:
    set val(sample_ID), file(bam_gatk_ra_rc) from sample_bam_delly2_inversions

    output:
    file "inversions.vcf"

    script:
    """
    $params.samtools_bin index ${bam_gatk_ra_rc}
    $params.delly2_bin call -t INV -g ${params.hg19_fa} -o inversions.bcf "${bam_gatk_ra_rc}"
    $params.delly2_bcftools_bin view inversions.bcf > inversions.vcf
    rm -f ${bam_gatk_ra_rc}.bai
    """
}

process delly2_translocations {
    tag { sample_ID }
    publishDir "${params.output_dir}/SNV-Delly2-translocations", mode: 'move', overwrite: true,
        saveAs: { filename -> "${sample_ID}_translocations.vcf" }

    input:
    set val(sample_ID), file(bam_gatk_ra_rc) from sample_bam_delly2_translocations

    output:
    file "translocations.vcf"

    script:
    """
    $params.samtools_bin index ${bam_gatk_ra_rc}
    $params.delly2_bin call -t BND -g ${params.hg19_fa} -o translocations.bcf "${bam_gatk_ra_rc}"
    $params.delly2_bcftools_bin view translocations.bcf > translocations.vcf
    rm -f ${bam_gatk_ra_rc}.bai
    """
}

process delly2_insertions {
    tag { sample_ID }
    publishDir "${params.output_dir}/SNV-Delly2-insertions", mode: 'move', overwrite: true,
        saveAs: { filename -> "${sample_ID}_insertions.vcf" }

    input:
    set val(sample_ID), file(bam_gatk_ra_rc) from sample_bam_delly2_insertions

    output:
    file "insertions.vcf"

    script:
    """
    $params.samtools_bin index ${bam_gatk_ra_rc}
    $params.delly2_bin call -t INS -g ${params.hg19_fa} -o insertions.bcf "${bam_gatk_ra_rc}"
    $params.delly2_bcftools_bin view insertions.bcf > insertions.vcf
    rm -f ${bam_gatk_ra_rc}.bai
    """
}


process  gatk_coverage_custom {
    tag { sample_ID }
    publishDir "${params.output_dir}/Coverage-GATK-Custom" , mode: 'move', overwrite: true

    input:
    file regions_bed name "regions.bed" from file(params.regions_bed)
    set val(sample_ID), file(bam_gatk_ra_rc) from sample_bam_gatk_coverage_custom

    output:
    file "${sample_ID}.sample_cumulative_coverage_counts"
    file "${sample_ID}.sample_cumulative_coverage_proportions"
    file "${sample_ID}.sample_interval_statistics"
    file "${sample_ID}.sample_interval_summary"
    file "${sample_ID}.sample_statistics"
    file "${sample_ID}.sample_summary"

    script:
    """
    $params.samtools_bin index ${bam_gatk_ra_rc}
    java -Xms16G -Xmx16G -jar ${params.gatk_bin} -T DepthOfCoverage \
    --logging_level ERROR \
    --downsampling_type NONE \
    --read_filter BadCigar \
    --reference_sequence ${params.hg19_fa} \
    --omitDepthOutputAtEachBase \
    -ct 10 -ct 50 -ct 100 -ct 200 -ct 300 -ct 400 -ct 500 \
    --intervals ${regions_bed} \
    --minBaseQuality 20 \
    --minMappingQuality 20 \
    --nBins 999 \
    --start 1 \
    --stop 1000 \
    --input_file ${bam_gatk_ra_rc} \
    --outputFormat csv \
    --out ${sample_ID}
    rm -f ${bam_gatk_ra_rc}.bai
    """
}
