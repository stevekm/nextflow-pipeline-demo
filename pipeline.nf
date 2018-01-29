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
                            def sample_tumor_ID = row."${params.tumor_header}"
                            def sample_tumor_bam = file("${params.bam_gatk_ra_rc_dir}/${sample_tumor_ID}.dd.ra.rc.bam")
                            def sample_tumor_bai = file("${params.bam_gatk_ra_rc_dir}/${sample_tumor_ID}.dd.ra.rc.bam.bai")
                            def sample_normal_ID = row."${params.normal_header}"
                            def sample_normal_bam = file("${params.bam_gatk_ra_rc_dir}/${sample_normal_ID}.dd.ra.rc.bam")
                            def sample_normal_bai = file("${params.bam_gatk_ra_rc_dir}/${sample_normal_ID}.dd.ra.rc.bam.bai")
                            def sample_ID = "${sample_tumor_ID}_${sample_normal_ID}" // sample comparison ID
                            return [ sample_ID, sample_tumor_ID, sample_tumor_bam, sample_tumor_bai, sample_normal_ID, sample_normal_bam, sample_normal_bai ]
                        }
                        .into {sample_pairs_demo;
                            sample_pairs_msi}

//
// read samples from variants summary sheet
//
Channel.fromPath( file(params.variants_sheet) )
                    .splitCsv(header: true)
                    .map{row ->
                        def sample_ID = row."${params.variant_sample_header}"
                        def sample_vcf = file("${params.variants_gatk_hc_dir}/${sample_ID}.vcf")
                        return [ sample_ID, sample_vcf ]
                    }
                    .set { sample_variants_deconstructSigs }



// Delly2 Channel; dont include the .bai, generate it instead
Channel.fromPath( file(params.variants_sheet) )
                    .splitCsv(header: true)
                    .map{row ->
                        def sample_ID = row."${params.variant_sample_header}"
                        def sample_bam = file("${params.bam_gatk_ra_rc_dir}/${sample_ID}.dd.ra.rc.bam")

                        println "-------------"
                        println "--- start csv row mapping ----"
                        println "row: ${row}"
                        println "sample_ID: ${sample_ID}"
                        println "sample_bam: ${sample_bam}"
                        println "--- end csv row mapping ----"
                        println "-------------"

                        return [ sample_ID, sample_bam ]
                    }
                    .into {
                        sample_bam_delly2_deletions;
                        sample_bam_delly2_duplications;
                        sample_bam_delly2_inversions;
                        sample_bam_delly2_translocations;
                        sample_bam_delly2_insertions;
                    }

// GATK channel; include .bai
Channel.fromPath( file(params.variants_sheet) )
                    .splitCsv(header: true)
                    .map{row ->
                        def sample_ID = row."${params.variant_sample_header}"
                        def sample_bam = file("${params.bam_gatk_ra_rc_dir}/${sample_ID}.dd.ra.rc.bam")
                        def sample_bai = file("${params.bam_gatk_ra_rc_dir}/${sample_ID}.dd.ra.rc.bam.bai")
                        return [ sample_ID, sample_bam, sample_bai ]
                    }
                    .into {
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
    set val(sample_ID), val(sample_tumor_ID), file(sample_tumor_bam), file(sample_tumor_bai), val(sample_normal_ID), file(sample_normal_bam), file(sample_normal_bai) from sample_pairs_demo

    exec:
    println "sample_ID: ${sample_ID}, sample_tumor_ID: ${sample_tumor_ID}, sample_tumor_bam: ${sample_tumor_bam}, sample_tumor_bai: ${sample_tumor_bai}, sample_normal_ID: ${sample_normal_ID}, sample_normal_bam: ${sample_normal_bam}, sample_normal_bai: ${sample_normal_bai}"

}


process msisensor {
    tag { sample_ID }
    clusterOptions '-pe threaded 1-4 -j y'
    publishDir "${params.output_dir}/MSI" , mode: 'move', overwrite: true, //"${params.output_dir}/${sample_ID}/MSI"
        saveAs: {filename ->
            if (filename == 'msisensor') "${sample_ID}.msisensor"
            else if (filename == 'msisensor_dis') "${sample_ID}.msisensor_dis"
            else if (filename == 'msisensor_germline') "${sample_ID}.msisensor_germline"
            else if (filename == 'msisensor_somatic') "${sample_ID}.msisensor_somatic"
            else null
        }

    input:
    set val(sample_ID), val(sample_tumor_ID), file(sample_tumor_bam), file(sample_tumor_bai), val(sample_normal_ID), file(sample_normal_bam), file(sample_normal_bai) from sample_pairs_msi
    file regions_bed from file(params.regions_bed) // name "regions.bed"
    // file microsatellites from file(params.microsatellites)

    output:
    file 'msisensor'
    file 'msisensor_dis'
    file 'msisensor_germline'
    file 'msisensor_somatic'

    // $params.samtools_bin index ${sample_normal_bam}
    // $params.samtools_bin index ${sample_tumor_bam}
    // rm -f ${sample_normal_bam}.bai
    // rm -f ${sample_tumor_bam}.bai
    script:
    """
    echo "USER:\${USER:-none}\tJOB_ID:\${JOB_ID:-none}\tJOB_NAME:\${JOB_NAME:-none}\tHOSTNAME:\${HOSTNAME:-none}\tTASK_ID:\${TASK_ID:-none}"
    $params.msisensor_bin msi -d ${params.microsatellites} -n $sample_normal_bam -t $sample_tumor_bam -e $regions_bed -o msisensor -l 1 -q 1 -b \${NSLOTS:-1}
    """
}


process check_samples_mapping {
    tag { sample_ID }
    executor "local"
    input:
    set val(sample_ID), file(sample_bam), file(sample_bai) from sample_bam_demo

    exec:
    println "sample_ID: ${sample_ID}, sample_bam: ${sample_bam}, sample_bai: ${sample_bai}"
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


process delly2_deletions {
    tag { sample_ID }
    publishDir "${params.output_dir}/SNV-Delly2-deletions", mode: 'move', overwrite: true,
        saveAs: { filename -> "${sample_ID}_deletions.vcf" }

    input:
    set val(sample_ID), file(sample_bam) from sample_bam_delly2_deletions

    output:
    file "deletions.vcf"

    script:
    """
    $params.samtools_bin index ${sample_bam}
    $params.delly2_bin call -t DEL -g ${params.hg19_fa} -o deletions.bcf "${sample_bam}"
    $params.delly2_bcftools_bin view deletions.bcf > deletions.vcf
    rm -f ${sample_bam}.bai
    """
}

process delly2_duplications {
    tag { sample_ID }
    publishDir "${params.output_dir}/SNV-Delly2-duplications", mode: 'move', overwrite: true,
        saveAs: { filename -> "${sample_ID}_duplications.vcf" }

    input:
    set val(sample_ID), file(sample_bam) from sample_bam_delly2_duplications

    output:
    file "duplications.vcf"

    script:
    """
    $params.samtools_bin index ${sample_bam}
    $params.delly2_bin call -t DUP -g ${params.hg19_fa} -o duplications.bcf "${sample_bam}"
    $params.delly2_bcftools_bin view duplications.bcf > duplications.vcf
    rm -f ${sample_bam}.bai
    """
}

process delly2_inversions {
    tag { sample_ID }
    publishDir "${params.output_dir}/SNV-Delly2-inversions", mode: 'move', overwrite: true,
        saveAs: { filename -> "${sample_ID}_inversions.vcf" }

    input:
    set val(sample_ID), file(sample_bam) from sample_bam_delly2_inversions

    output:
    file "inversions.vcf"

    script:
    """
    $params.samtools_bin index ${sample_bam}
    $params.delly2_bin call -t INV -g ${params.hg19_fa} -o inversions.bcf "${sample_bam}"
    $params.delly2_bcftools_bin view inversions.bcf > inversions.vcf
    rm -f ${sample_bam}.bai
    """
}

process delly2_translocations {
    tag { sample_ID }
    publishDir "${params.output_dir}/SNV-Delly2-translocations", mode: 'move', overwrite: true,
        saveAs: { filename -> "${sample_ID}_translocations.vcf" }

    input:
    set val(sample_ID), file(sample_bam) from sample_bam_delly2_translocations

    output:
    file "translocations.vcf"

    script:
    """
    $params.samtools_bin index ${sample_bam}
    $params.delly2_bin call -t BND -g ${params.hg19_fa} -o translocations.bcf "${sample_bam}"
    $params.delly2_bcftools_bin view translocations.bcf > translocations.vcf
    rm -f ${sample_bam}.bai
    """
}

process delly2_insertions {
    tag { sample_ID }
    publishDir "${params.output_dir}/SNV-Delly2-insertions", mode: 'move', overwrite: true,
        saveAs: { filename -> "${sample_ID}_insertions.vcf" }

    input:
    set val(sample_ID), file(sample_bam) from sample_bam_delly2_insertions

    output:
    file "insertions.vcf"

    script:
    """
    $params.samtools_bin index ${sample_bam}
    $params.delly2_bin call -t INS -g ${params.hg19_fa} -o insertions.bcf "${sample_bam}"
    $params.delly2_bcftools_bin view insertions.bcf > insertions.vcf
    rm -f ${sample_bam}.bai
    """
}


process  gatk_coverage_custom {
    tag { sample_ID }
    publishDir "${params.output_dir}/Coverage-GATK-Custom" , mode: 'move', overwrite: true

    input:
    file regions_bed from file(params.regions_bed) // name "regions.bed"
    set val(sample_ID), file(sample_bam), file(sample_bai) from sample_bam_gatk_coverage_custom

    output:
    file "${sample_ID}.sample_cumulative_coverage_counts"
    file "${sample_ID}.sample_cumulative_coverage_proportions"
    file "${sample_ID}.sample_interval_statistics"
    file "${sample_ID}.sample_interval_summary"
    file "${sample_ID}.sample_statistics"
    file "${sample_ID}.sample_summary"

    script:
    """
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
    --input_file ${sample_bam} \
    --outputFormat csv \
    --out ${sample_ID}
    """
}
