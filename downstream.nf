// params.input_dir = "sns-dir"
// params.regions_bed = "targets.bed"
// params.probes_bed = "probes.bed"
//
// params.bam_bwa_dir = "${params.input_dir}/BAM-BWA"
// params.bam_gatk_ra_rc_dir = "${params.input_dir}/BAM-GATK-RA-RC"
//
// params.pairs_sheet = "${params.input_dir}/samples.pairs.csv"
// params.variants_sheet = "${params.input_dir}/summary.VCF-GATK-HC-annot.csv"
// params.variants_gatk_hc_dir = "${params.input_dir}/VCF-GATK-HC"
//
// params.variants_gatk_hc_annot_dir = "${params.input_dir}/VCF-GATK-HC-annot"
// // /ifs/data/molecpathlab/PNET_GYN/sns_WES/VCF-GATK-HC-annot/256.combined.txt
//
// params.qc_coverage_dir = "${params.input_dir}/QC-coverage"
// // sns-dir/QC-coverage/*.sample_interval_summary

// get full path to the script
// import java.io.File;
r_util_script = new File(params.r_util_script).getCanonicalPath()


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
                            return [ sample_ID, sample_tumor_ID, sample_tumor_bam,
                                    // sample_tumor_bai,
                                    sample_normal_ID, sample_normal_bam
                                    // sample_normal_bai
                                    ]
                        }
                        .filter{ sample_ID, sample_tumor_ID, sample_tumor_bam, sample_normal_ID, sample_normal_bam ->
                            sample_tumor_ID != params.pairs_sheet_nomatch_value && sample_normal_ID != params.pairs_sheet_nomatch_value
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
                        sample_bam_gatk_coverage_custom;
                        sample_bam_demo
                    }

// GATK channel; include .bai
// Channel.fromPath( file(params.variants_sheet) )
//                     .splitCsv(header: true)
//                     .map{row ->
//                         def sample_ID = row."${params.variant_sample_header}"
//                         def sample_bam = file("${params.bam_gatk_ra_rc_dir}/${sample_ID}.dd.ra.rc.bam")
//                         def sample_bai = file("${params.bam_gatk_ra_rc_dir}/${sample_ID}.dd.ra.rc.bam.bai")
//                         return [ sample_ID, sample_bam, sample_bai ]
//                     }
//                     .into {
//                         sample_bam_gatk_coverage_custom;
//                         sample_bam_demo
//                     }

// VAF plot channel
Channel.fromPath( file(params.variants_sheet) )
                    .splitCsv(header: true)
                    .map{row ->
                        def sample_ID = row."${params.variant_sample_header}"
                        def sample_vcf_annot = file("${params.variants_gatk_hc_annot_dir}/${sample_ID}.combined.txt")
                        return [ sample_ID, sample_vcf_annot ]
                    }
                    .into {
                        sample_vcf_annot;
                        sample_vcf_annot2
                    }


// files to send in email
email_files = Channel.create()


//
// pipeline steps
//


// DEBUGGING STEPS
process match_samples {
    tag { sample_ID }
    executor "local"
    input:
    set val(sample_ID), val(sample_tumor_ID), file(sample_tumor_bam), val(sample_normal_ID), file(sample_normal_bam) from sample_pairs_demo
    // set val(sample_ID), val(sample_tumor_ID), file(sample_tumor_bam), file(sample_tumor_bai), val(sample_normal_ID), file(sample_normal_bam), file(sample_normal_bai) from sample_pairs_demo

    exec:
    println "sample_ID: ${sample_ID}, sample_tumor_ID: ${sample_tumor_ID}, sample_tumor_bam: ${sample_tumor_bam}, sample_normal_ID: ${sample_normal_ID}, sample_normal_bam: ${sample_normal_bam}"
    // println "sample_ID: ${sample_ID}, sample_tumor_ID: ${sample_tumor_ID}, sample_tumor_bam: ${sample_tumor_bam}, sample_tumor_bai: ${sample_tumor_bai}, sample_normal_ID: ${sample_normal_ID}, sample_normal_bam: ${sample_normal_bam}, sample_normal_bai: ${sample_normal_bai}"

}


process check_samples_mapping {
    tag { sample_ID }
    executor "local"
    input:
    set val(sample_ID), file(sample_bam) from sample_bam_demo

    exec:
    println "sample_ID: ${sample_ID}, sample_bam: ${sample_bam}"
}


// NEW PIPELINE STEPS
// process deconstructSigs_signatures {
//     tag { sample_ID }
//     // executor "local"
//     validExitStatus 0,11 // allow '11' failure triggered by no variants
//     errorStrategy 'ignore'
//     publishDir "${params.output_dir}/Genomic_Signatures", mode: 'copy', overwrite: true,
//         saveAs: {filename ->
//             if (filename == 'sample_signatures.Rds') "${sample_ID}_signatures.Rds"
//             else if (filename == 'sample_signatures.pdf') "${sample_ID}_signatures.pdf"
//             else if (filename == 'sample_signatures_pie.pdf') "${sample_ID}_signatures_pie.pdf"
//             else filename
//         }
//     input:
//     set val(sample_ID), file(sample_vcf) from sample_variants_deconstructSigs
//
//     output:
//     file "${sample_ID}_signatures.Rds"
//     file "${sample_ID}_signatures.pdf" into signatures_plots
//     file "${sample_ID}_signatures_pie.pdf" into signatures_pie_plots
//
//     script:
//     """
//     pwd
//     $params.deconstructSigs_make_signatures_script "${sample_ID}" "$sample_vcf"
//     """
// }
//
//
// process vaf_distribution_plot {
//     tag { sample_ID }
//     // executor "local"
//     publishDir "${params.output_dir}/VAF-Distribution", mode: 'copy', overwrite: true
//
//     input:
//     set val(sample_ID), file(sample_vcf_annot) from sample_vcf_annot
//
//     output:
//     file "${sample_ID}_vaf_dist.pdf" into vaf_distribution_plots
//
//     script:
//     """
//     $params.vaf_distribution_plot_script "${sample_ID}" "${sample_vcf_annot}"
//     """
//
// }


//
// process msisensor {
//     tag { sample_ID }
//     module 'samtools/1.3'
//     clusterOptions '-pe threaded 1-4 -j y -l mem_free=40G'
//     publishDir "${params.output_dir}/MSI", mode: 'copy', overwrite: true,
//         saveAs: {filename ->
//             if (filename == 'msisensor') "${sample_ID}.msisensor"
//             else if (filename == 'msisensor_dis') "${sample_ID}.msisensor_dis"
//             else if (filename == 'msisensor_germline') "${sample_ID}.msisensor_germline"
//             else if (filename == 'msisensor_somatic') "${sample_ID}.msisensor_somatic"
//             else null
//         }
//
//     input:
//     set val(sample_ID), val(sample_tumor_ID), file(sample_tumor_bam), val(sample_normal_ID), file(sample_normal_bam) from sample_pairs_msi
//     file regions_bed from file(params.regions_bed) // name "regions.bed"
//     // set val(sample_ID), val(sample_tumor_ID), file(sample_tumor_bam), file(sample_tumor_bai), val(sample_normal_ID), file(sample_normal_bam), file(sample_normal_bai) from sample_pairs_msi
//     // file microsatellites from file(params.microsatellites)
//
//     output:
//     file 'msisensor'
//     file 'msisensor_dis'
//     file 'msisensor_germline'
//     file 'msisensor_somatic'
//
//     script:
//     // $params.subset_msisensor_microsatellite_list_script *.bam -m $params.microsatellites -o microsatellites_subset.txt
//     // $params.subset_bam_bed_script *.bam -b $regions_bed -o targets_subset.bed
//     // rm -f "${sample_normal_bam}.bai" "${sample_tumor_bam}.bai" targets_subset.bed microsatellites_subset.txt
//     """
//     echo "USER:\${USER:-none}\tJOB_ID:\${JOB_ID:-none}\tJOB_NAME:\${JOB_NAME:-none}\tHOSTNAME:\${HOSTNAME:-none}\tTASK_ID:\${TASK_ID:-none}"
//     samtools index "$sample_tumor_bam"
//     samtools index "$sample_normal_bam"
//     $params.msisensor_bin msi -d $params.microsatellites -n $sample_normal_bam -t $sample_tumor_bam -e $regions_bed -o msisensor -l 1 -q 1 -b \${NSLOTS:-1}
//     rm -f "${sample_normal_bam}.bai" "${sample_tumor_bam}.bai"
//     """
// }


// OLD PIPELINE STEPS
// // DELLY2 SNV STEPS
// process delly2_deletions {
//     tag { sample_ID }
//     publishDir "${params.output_dir}/SNV-Delly2-deletions", mode: 'copy', overwrite: true,
//         saveAs: { filename -> "${sample_ID}_deletions.vcf" }
//
//     input:
//     set val(sample_ID), file(sample_bam) from sample_bam_delly2_deletions
//
//     output:
//     file "deletions.vcf"
//
//     script:
//     """
//     $params.samtools_bin index ${sample_bam}
//     $params.delly2_bin call -t DEL -g ${params.hg19_fa} -o deletions.bcf "${sample_bam}"
//     $params.delly2_bcftools_bin view deletions.bcf > deletions.vcf
//     rm -f ${sample_bam}.bai
//     """
// }
//
// process delly2_duplications {
//     tag { sample_ID }
//     publishDir "${params.output_dir}/SNV-Delly2-duplications", mode: 'copy', overwrite: true,
//         saveAs: { filename -> "${sample_ID}_duplications.vcf" }
//
//     input:
//     set val(sample_ID), file(sample_bam) from sample_bam_delly2_duplications
//
//     output:
//     file "duplications.vcf"
//
//     script:
//     """
//     $params.samtools_bin index ${sample_bam}
//     $params.delly2_bin call -t DUP -g ${params.hg19_fa} -o duplications.bcf "${sample_bam}"
//     $params.delly2_bcftools_bin view duplications.bcf > duplications.vcf
//     rm -f ${sample_bam}.bai
//     """
// }
//
// process delly2_inversions {
//     tag { sample_ID }
//     publishDir "${params.output_dir}/SNV-Delly2-inversions", mode: 'copy', overwrite: true,
//         saveAs: { filename -> "${sample_ID}_inversions.vcf" }
//
//     input:
//     set val(sample_ID), file(sample_bam) from sample_bam_delly2_inversions
//
//     output:
//     file "inversions.vcf"
//
//     script:
//     """
//     $params.samtools_bin index ${sample_bam}
//     $params.delly2_bin call -t INV -g ${params.hg19_fa} -o inversions.bcf "${sample_bam}"
//     $params.delly2_bcftools_bin view inversions.bcf > inversions.vcf
//     rm -f ${sample_bam}.bai
//     """
// }
//
// process delly2_translocations {
//     tag { sample_ID }
//     publishDir "${params.output_dir}/SNV-Delly2-translocations", mode: 'copy', overwrite: true,
//         saveAs: { filename -> "${sample_ID}_translocations.vcf" }
//
//     input:
//     set val(sample_ID), file(sample_bam) from sample_bam_delly2_translocations
//
//     output:
//     file "translocations.vcf"
//
//     script:
//     """
//     $params.samtools_bin index ${sample_bam}
//     $params.delly2_bin call -t BND -g ${params.hg19_fa} -o translocations.bcf "${sample_bam}"
//     $params.delly2_bcftools_bin view translocations.bcf > translocations.vcf
//     rm -f ${sample_bam}.bai
//     """
// }
//
// process delly2_insertions {
//     tag { sample_ID }
//     publishDir "${params.output_dir}/SNV-Delly2-insertions", mode: 'copy', overwrite: true,
//         saveAs: { filename -> "${sample_ID}_insertions.vcf" }
//
//     input:
//     set val(sample_ID), file(sample_bam) from sample_bam_delly2_insertions
//
//     output:
//     file "insertions.vcf"
//
//     script:
//     """
//     $params.samtools_bin index ${sample_bam}
//     $params.delly2_bin call -t INS -g ${params.hg19_fa} -o insertions.bcf "${sample_bam}"
//     $params.delly2_bcftools_bin view insertions.bcf > insertions.vcf
//     rm -f ${sample_bam}.bai
//     """
// }



process  gatk_coverage_custom {
    tag { sample_ID }
    module 'samtools/1.3'
    publishDir "${params.output_dir}/Coverage-GATK-Custom", mode: 'copy', overwrite: true

    input:
    file regions_bed from file(params.regions_bed) // name "regions.bed"
    set val(sample_ID), file(sample_bam) from sample_bam_gatk_coverage_custom

    output:
    file "${sample_ID}.sample_cumulative_coverage_counts"
    file "${sample_ID}.sample_cumulative_coverage_proportions"
    file "${sample_ID}.sample_interval_statistics"
    file "${sample_ID}.sample_interval_summary" into sample_interval_summary
    file "${sample_ID}.sample_statistics"
    file "${sample_ID}.sample_summary"

    script:
    """
    $params.samtools_bin index ${sample_bam}

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

    rm -f ${sample_bam}.bai
    """
}



// CUSTOM PLOTTING & SCRIPT STEPS


process summary_GATK_intervals {
    echo true
    executor "local"
    publishDir "${params.output_dir}/Coverage-GATK-Interval-Summary", mode: 'copy', overwrite: true

    input:
    file "*" from sample_interval_summary.toList()

    output:
    file "average_coverage_per_sample.tsv"
    file "average_coverage_per_region.tsv"
    file "regions_coverage_below_50.bed"
    file "regions_with_coverage_0.bed"

    script:
    """
    set -x
    echo "running summary_GATK_intervals task..."
    echo "pwd is: \$(pwd)"
    $params.gatk_avg_coverages_script *.sample_interval_summary && echo "done script on process: \$\$, exit status: \$?" && exit 0
    """
}


process merge_VAF_plots {
    executor "local"
    publishDir "${params.output_dir}/VAF-Distribution_Summary", mode: 'copy', overwrite: true

    input:
    file '*' from vaf_distribution_plots.toList()

    output:
    file "vaf_distributions.pdf" into email_VAF_plots

    script:
    """
    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=vaf_distributions.pdf *
    """
}

process merge_signatures_plots {
    validExitStatus 0,11 // allow '11' failure triggered by no variants
    errorStrategy 'ignore'
    executor "local"
    publishDir "${params.output_dir}/Genomic_Signatures_Summary", mode: 'copy', overwrite: true

    input:
    file '*' from signatures_plots.toList()

    output:
    file "genomic_signatures.pdf" into email_signatures_plots

    script:
    """
    if [ "\$(ls -1 * | wc -l)" -gt 0 ]; then
        gs -dBATCH -dNOPAUSE -q -dAutoRotatePages=/None -sDEVICE=pdfwrite -sOutputFile=genomic_signatures.pdf *
    else
        exit 11
    fi
    """
}


process merge_signatures_pie_plots {
    validExitStatus 0,11 // allow '11' failure triggered by no variants
    errorStrategy 'ignore'
    executor "local"
    publishDir "${params.output_dir}/Genomic_Signatures_Summary", mode: 'copy', overwrite: true

    input:
    file '*' from signatures_pie_plots.toList()

    output:
    file "genomic_signatures_pie.pdf" into email_signatures_pie_plots

    script:
    """
    if [ "\$(ls -1 * | wc -l)" -gt 0 ]; then
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=genomic_signatures_pie.pdf *
    else
        exit 11
    fi
    """
}

//
//
// process email {
//     echo true
//     executor "local"
//
//     input:
//     file '*' from email_files.mix(email_VAF_plots,
//                                     email_signatures_pie_plots,
//                                     email_signatures_plots
//                                     ).toList()
//
//     script:
//     """
//     echo *
//     """
// }
