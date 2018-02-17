// per-sample exome analysis
params.output_dir = "output-exome"

// read samples from fastq raw sheet, group R1 and R2 files per sample
Channel.fromPath( file(params.fastq_raw_sheet) )
        .splitCsv()
        .map { row ->
            def sample_ID = row[0]
            def read1 = file(row[1])
            def read2 = file(row[2])
            return [ sample_ID, read1, read2 ]
        }
        .groupTuple()
        .into { sample_fastq_r1r2; sample_fastq_r1r2_2 }

// target regions .bed file
Channel.fromPath( file(params.targets_bed) )
        .into { targets_bed;
                targets_bed2;
                targets_bed3;
                targets_bed4;
                targets_bed5;
                targets_bed6;
                targets_bed7;
                targets_bed8 }

// reference hg19 fasta file
Channel.fromPath( file(params.hg19_fa) )
        .into { ref_fasta;
                ref_fasta2;
                ref_fasta3;
                ref_fasta4;
                ref_fasta5;
                ref_fasta6;
                ref_fasta7;
                ref_fasta8;
                ref_fasta9 }
Channel.fromPath( file(params.hg19_fai) )
        .into { ref_fai;
                ref_fai2;
                ref_fai3;
                ref_fai4;
                ref_fai5;
                ref_fai6;
                ref_fai7;
                ref_fai8;
                ref_fai9 }
Channel.fromPath( file(params.hg19_dict) )
        .into { ref_dict;
                ref_dict2;
                ref_dict3;
                ref_dict4;
                ref_dict5;
                ref_dict6;
                ref_dict7;
                ref_dict8;
                ref_dict9 }

Channel.fromPath( file(params.hg19_chrom_sizes) )
        .set { ref_chrom_sizes }

//
//
// DEBUGGING
//
//


process fastq_pairs_print {
    // prints the fastq file R1 and R2 sets
    executor "local"
    echo true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"

    input:
    set val(sample_ID), file(fastq_r1: "*"), file(fastq_r2: "*") from sample_fastq_r1r2_2

    script:
    """
    echo "----------------"
    echo "${sample_ID} - \$(pwd)"
    echo "R1 files: ${fastq_r1}"
    echo "R2 files: ${fastq_r2}"
    echo "----------------"
    """
}

//
//
// EXOME SEQUENCING ANALYSIS PIPELINE
//
//
//
// process fastq_merge {
//     // merge the R1 and R2 fastq files into a single fastq each
//     tag { "${sample_ID}" }
//     beforeScript "${params.beforeScript_str}"
//     afterScript "${params.afterScript_str}"
//     publishDir "${params.wes_output_dir}/fastq-merge", mode: 'copy', overwrite: true
//
//     input:
//     set val(sample_ID), file(fastq_r1: "*"), file(fastq_r2: "*") from sample_fastq_r1r2
//
//     output:
//     set val(sample_ID), file("${sample_ID}_R1.fastq.gz"), file("${sample_ID}_R2.fastq.gz") into samples_fastq_merged
//
//     script:
//     """
//     cat ${fastq_r1} > "${sample_ID}_R1.fastq.gz"
//     cat ${fastq_r2} > "${sample_ID}_R2.fastq.gz"
//     """
// }

//
// process trimmomatic {
//     // Illumina read trimming
//     // http://www.usadellab.org/cms/?page=trimmomatic
//     tag { "${sample_ID}-${read1}-${read2}" }
//     publishDir "${params.wes_output_dir}/fastq-trim", mode: 'copy', overwrite: true
//     clusterOptions '-pe threaded 1-8'
//     beforeScript "${params.beforeScript_str}"
//     afterScript "${params.afterScript_str}"
//
//     input:
//     set val(sample_ID), file(read1), file(read2) from samples_fastq_merged
//
//     output:
//     set val(sample_ID), file("${sample_ID}_R1.trim.fastq.gz"), file("${sample_ID}_R1.unpaired.fastq.gz"), file("${sample_ID}_R2.trim.fastq.gz"), file("${sample_ID}_R2.unpaired.fastq.gz") into samples_fastq_trimmed
//
//     script:
//     """
//     # echo "${sample_ID} ${read1} ${read2} \${NSLOTS:-1}"
//     java -Xms16G -Xmx16G -jar ${params.trimmomatic_jar} PE -threads \${NSLOTS:-1} \
//     "${read1}" "${read2}" \
//     "${sample_ID}_R1.trim.fastq.gz" "${sample_ID}_R1.unpaired.fastq.gz" \
//     "${sample_ID}_R2.trim.fastq.gz" "${sample_ID}_R2.unpaired.fastq.gz" \
//     ILLUMINACLIP:${params.trimmomatic_contaminant_fa}:2:30:10:1:true TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:35
//     """
// }
//
// process bwa_mem {
//     // first pass alignment with BWA
//     tag { "${sample_ID}" }
//     clusterOptions '-pe threaded 1-8'
//     beforeScript "${params.beforeScript_str}"
//     afterScript "${params.afterScript_str}"
//     module 'bwa/0.7.17'
//
//     input:
//     set val(sample_ID), file(fastq_R1_trim), file(fastq_R1_trim_unpaired), file(fastq_R2_trim), file(fastq_R2_trim_unpaired) from samples_fastq_trimmed
//
//     output:
//     set val(sample_ID), file("${sample_ID}.sam") into samples_bwa_sam
//
//     script:
//     """
//     bwa mem -M -v 1 -t \${NSLOTS:-1} -R '@RG\\tID:${sample_ID}\\tSM:${sample_ID}\\tLB:${sample_ID}\\tPL:ILLUMINA' "${params.bwa_hg19_ref_fa}" "${fastq_R1_trim}" "${fastq_R2_trim}" -o "${sample_ID}.sam"
//     """
// }
//
// process sambamba_view_sort {
//     tag { "${sample_ID}" }
//     clusterOptions '-pe threaded 1-8'
//     beforeScript "${params.beforeScript_str}"
//     afterScript "${params.afterScript_str}"
//
//     input:
//     set val(sample_ID), file(sample_sam) from samples_bwa_sam
//
//     output:
//     set val(sample_ID), file("${sample_ID}.bam") into samples_bam, samples_bam2
//
//     script:
//     """
//     "${params.sambamba_bin}" view --sam-input --nthreads=\${NSLOTS:-1} --filter='mapping_quality>=10' --format=bam --compression-level=0 "${sample_sam}" | \
//     "${params.sambamba_bin}" sort --nthreads=\${NSLOTS:-1} --memory-limit=16GB --out="${sample_ID}.bam" /dev/stdin
//     """
// }
//
// process sambamba_flagstat {
//     tag { "${sample_ID}" }
//     publishDir "${params.wes_output_dir}/sambamba-flagstat", mode: 'copy', overwrite: true
//     beforeScript "${params.beforeScript_str}"
//     afterScript "${params.afterScript_str}"
//
//     input:
//     set val(sample_ID), file(sample_bam) from samples_bam2
//
//     output:
//     file "${sample_ID}.flagstat.txt"
//
//     script:
//     """
//     "${params.sambamba_bin}" flagstat "${sample_bam}" > "${sample_ID}.flagstat.txt"
//     """
// }
//
// process sambamba_dedup {
//     tag { "${sample_ID}" }
//     publishDir "${params.wes_output_dir}/bam-bwa", mode: 'copy', overwrite: true
//     clusterOptions '-pe threaded 1-8'
//     beforeScript "${params.beforeScript_str}"
//     afterScript "${params.afterScript_str}"
//
//     input:
//     set val(sample_ID), file(sample_bam) from samples_bam
//
//     output:
//     set val(sample_ID), file("${sample_ID}.dd.bam") into samples_dd_bam, samples_dd_bam2, samples_dd_bam3, samples_dd_bam4, samples_dd_bam5, samples_dd_bam6, samples_dd_bam7
//
//     script:
//     """
//     "${params.sambamba_bin}" markdup --remove-duplicates --nthreads \${NSLOTS:-1} --hash-table-size 525000 --overflow-list-size 525000 "${sample_bam}" "${sample_ID}.dd.bam"
//     """
// }
//
// process sambamba_dedup_flagstat {
//     tag { "${sample_ID}" }
//     // executor "local"
//     publishDir "${params.wes_output_dir}/sambamba-flagstat-dd", mode: 'copy', overwrite: true
//     beforeScript "${params.beforeScript_str}"
//     afterScript "${params.afterScript_str}"
//
//     input:
//     set val(sample_ID), file(sample_bam) from samples_dd_bam2
//
//     output:
//     file "${sample_ID}.dd.flagstat.txt"
//
//     script:
//     """
//     "${params.sambamba_bin}" flagstat "${sample_bam}" > "${sample_ID}.dd.flagstat.txt"
//     """
//
// }
//
// process qc_target_reads_gatk_genome {
//     tag { "${sample_ID}" }
//     // executor "local"
//     publishDir "${params.wes_output_dir}/qc-target-reads", mode: 'copy', overwrite: true
//     beforeScript "${params.beforeScript_str}"
//     afterScript "${params.afterScript_str}"
//     clusterOptions '-pe threaded 1-8'
//
//     input:
//     set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict) from samples_dd_bam.combine(ref_fasta).combine(ref_fai).combine(ref_dict)
//
//     output:
//     file "${sample_ID}.genome.sample_statistics"
//     file "${sample_ID}.genome.sample_summary"
//
//     script:
//     """
//     java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage -dt NONE -rf BadCigar -nt \${NSLOTS:-1} --logging_level ERROR --omitIntervalStatistics --omitLocusTable --omitDepthOutputAtEachBase -ct 10 -ct 100 -mbq 20 -mmq 20 --reference_sequence "${ref_fasta}" --input_file "${sample_bam}" --outputFormat csv --out "${sample_ID}.genome"
//     """
//
// }
//
// process qc_target_reads_gatk_pad500 {
//     tag { "${sample_ID}" }
//     // executor "local"
//     publishDir "${params.wes_output_dir}/qc-target-reads", mode: 'copy', overwrite: true
//     beforeScript "${params.beforeScript_str}"
//     afterScript "${params.afterScript_str}"
//     clusterOptions '-pe threaded 1-8'
//
//     input:
//     set val(sample_ID), file(sample_bam), file(targets_bed_file), file(ref_fasta), file(ref_fai), file(ref_dict) from samples_dd_bam3.combine(targets_bed).combine(ref_fasta2).combine(ref_fai2).combine(ref_dict2)
//
//     output:
//     file "${sample_ID}.pad500.sample_statistics"
//     file "${sample_ID}.pad500.sample_summary"
//
//     script:
//     """
//     java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage -dt NONE -rf BadCigar -nt \${NSLOTS:-1} --logging_level ERROR --omitIntervalStatistics --omitLocusTable --omitDepthOutputAtEachBase -ct 10 -ct 100 -mbq 20 -mmq 20 --reference_sequence "${ref_fasta}" --intervals "${targets_bed_file}" --interval_padding 500 --input_file "${sample_bam}" --outputFormat csv --out "${sample_ID}.pad500"
//     """
// }
//
// process qc_target_reads_gatk_pad100 {
//     tag { "${sample_ID}" }
//     publishDir "${params.wes_output_dir}/qc-target-reads", mode: 'copy', overwrite: true
//     beforeScript "${params.beforeScript_str}"
//     afterScript "${params.afterScript_str}"
//     clusterOptions '-pe threaded 1-8'
//
//     input:
//     set val(sample_ID), file(sample_bam), file(targets_bed_file), file(ref_fasta), file(ref_fai), file(ref_dict) from samples_dd_bam4.combine(targets_bed2).combine(ref_fasta3).combine(ref_fai3).combine(ref_dict3)
//
//     output:
//     file "${sample_ID}.pad100.sample_statistics"
//     file "${sample_ID}.pad100.sample_summary"
//
//     script:
//     """
//     java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage -dt NONE -rf BadCigar -nt \${NSLOTS:-1} --logging_level ERROR --omitIntervalStatistics --omitLocusTable --omitDepthOutputAtEachBase -ct 10 -ct 100 -mbq 20 -mmq 20 --reference_sequence "${ref_fasta}" --intervals "${targets_bed_file}" --interval_padding 100 --input_file "${sample_bam}" --outputFormat csv --out "${sample_ID}.pad100"
//     """
// }
//
// process qc_target_reads_gatk_bed {
//     tag { "${sample_ID}" }
//     publishDir "${params.wes_output_dir}/qc-target-reads", mode: 'copy', overwrite: true
//     beforeScript "${params.beforeScript_str}"
//     afterScript "${params.afterScript_str}"
//     clusterOptions '-pe threaded 1-8'
//
//     input:
//     set val(sample_ID), file(sample_bam), file(targets_bed_file), file(ref_fasta), file(ref_fai), file(ref_dict) from samples_dd_bam5.combine(targets_bed3).combine(ref_fasta4).combine(ref_fai4).combine(ref_dict4)
//
//     output:
//     file "${sample_ID}.bed.sample_statistics"
//     file "${sample_ID}.bed.sample_summary"
//
//     script:
//     """
//     java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage -dt NONE -rf BadCigar -nt \${NSLOTS:-1} --logging_level ERROR --omitIntervalStatistics --omitLocusTable --omitDepthOutputAtEachBase -ct 10 -ct 100 -mbq 20 -mmq 20 --reference_sequence "${ref_fasta}" --intervals "${targets_bed_file}" --input_file "${sample_bam}" --outputFormat csv --out "${sample_ID}.bed"
//     """
// }

process bam_ra_rc_gatk {
    // re-alignment and recalibration with GATK
    // https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
    // https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php
    tag { "${sample_ID}" }
    publishDir "${params.wes_output_dir}/${params.bam_ra_rc_gatk_dir}", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'

    input:
    set val(sample_ID), file(sample_bam), file(targets_bed_file), file(ref_fasta), file(ref_fai), file(ref_dict) from samples_dd_bam6.combine(targets_bed4).combine(ref_fasta5).combine(ref_fai5).combine(ref_dict5)

    output:
    set val(sample_ID), file("${sample_ID}.dd.ra.rc.bam") into samples_dd_ra_rc_bam, samples_dd_ra_rc_bam2, samples_dd_ra_rc_bam3, samples_dd_ra_rc_bam4
    file "${sample_ID}.intervals"
    file "${sample_ID}.table1.txt"
    file "${sample_ID}.table2.txt"
    file "${sample_ID}.csv"
    file "${sample_ID}.pdf"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T RealignerTargetCreator \
    -dt NONE \
    --logging_level ERROR \
    -nt \${NSLOTS:-1} \
    --reference_sequence "${ref_fasta}" \
    -known "${params.gatk_1000G_phase1_indels_hg19_vcf}" \
    -known "${params.mills_and_1000G_gold_standard_indels_hg19_vcf}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${sample_bam}" \
    --out "${sample_ID}.intervals"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T IndelRealigner \
    -dt NONE \
    --logging_level ERROR \
    --reference_sequence "${ref_fasta}" \
    --maxReadsForRealignment 50000 \
    -known "${params.gatk_1000G_phase1_indels_hg19_vcf}" \
    -known "${params.mills_and_1000G_gold_standard_indels_hg19_vcf}" \
    -targetIntervals "${sample_ID}.intervals" \
    --input_file "${sample_bam}" \
    --out "${sample_ID}.dd.ra.bam"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T BaseRecalibrator \
    --logging_level ERROR \
    -nct \${NSLOTS:-1} \
    -rf BadCigar \
    --reference_sequence "${ref_fasta}" \
    -knownSites "${params.gatk_1000G_phase1_indels_hg19_vcf}" \
    -knownSites "${params.mills_and_1000G_gold_standard_indels_hg19_vcf}" \
    -knownSites "${params.dbsnp_138_hg19_vcf}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${sample_ID}.dd.ra.bam" \
    --out "${sample_ID}.table1.txt"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T BaseRecalibrator \
    --logging_level ERROR \
    -nct \${NSLOTS:-1} \
    -rf BadCigar \
    --reference_sequence "${ref_fasta}" \
    -knownSites "${params.gatk_1000G_phase1_indels_hg19_vcf}" \
    -knownSites "${params.mills_and_1000G_gold_standard_indels_hg19_vcf}" \
    -knownSites "${params.dbsnp_138_hg19_vcf}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${sample_ID}.dd.ra.bam" \
    -BQSR "${sample_ID}.table1.txt" \
    --out "${sample_ID}.table2.txt"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T AnalyzeCovariates \
    --logging_level ERROR \
    --reference_sequence "${ref_fasta}" \
    -before "${sample_ID}.table1.txt" \
    -after "${sample_ID}.table2.txt" \
    -csv "${sample_ID}.csv" \
    -plots "${sample_ID}.pdf"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T PrintReads \
    --logging_level ERROR \
    -nct \${NSLOTS:-1} \
    -rf BadCigar \
    --reference_sequence "${ref_fasta}" \
    -BQSR "${sample_ID}.table1.txt" \
    --input_file "${sample_ID}.dd.ra.bam" \
    --out "${sample_ID}.dd.ra.rc.bam"
    """

}

process qc_coverage_gatk {
    tag { "${sample_ID}" }
    publishDir "${params.wes_output_dir}/qc_coverage_gatk", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'

    input:
    set val(sample_ID), file(sample_bam), file(targets_bed_file), file(ref_fasta), file(ref_fai), file(ref_dict) from samples_dd_ra_rc_bam.combine(targets_bed5).combine(ref_fasta6).combine(ref_fai6).combine(ref_dict6)

    output:
    file "${sample_ID}.sample_summary"
    file "${sample_ID}.sample_statistics"
    file "${sample_ID}.sample_interval_summary"
    file "${sample_ID}.sample_interval_statistics"
    file "${sample_ID}.sample_cumulative_coverage_proportions"
    file "${sample_ID}.sample_cumulative_coverage_counts"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage \
    -dt NONE \
    --logging_level ERROR \
    -rf BadCigar \
    --reference_sequence "${ref_fasta}" \
    --intervals "${targets_bed_file}" \
    --omitDepthOutputAtEachBase \
    -ct 10 -ct 50 -ct 100 -ct 500 \
    -mbq 20 -mmq 20 --nBins 999 \
    --start 1 --stop 1000 \
    --input_file "${sample_bam}" \
    --outputFormat csv \
    --out "${sample_ID}"
    """
}

process pad_bed {
    publishDir "${params.wes_output_dir}/targets", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    module 'bedtools/2.26.0'

    input:
    set file(targets_bed_file), file(ref_chrom_sizes) from targets_bed6.combine(ref_chrom_sizes)

    output:
    file("targets.pad10.bed") into targets_pad_bed

    script:
    """
    cat "${targets_bed_file}" | LC_ALL=C sort -k1,1 -k2,2n | bedtools slop -g "${ref_chrom_sizes}" -b 10 | bedtools merge -d 5 > targets.pad10.bed
    """
}

process lofreq {
    tag { "${sample_ID}" }
    publishDir "${params.wes_output_dir}/vcf_lofreq", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'
    module 'samtools/1.3'

    input:
    set val(sample_ID), file(sample_bam), file(targets_bed_file), file(ref_fasta), file(ref_fai) from samples_dd_ra_rc_bam2.combine(targets_pad_bed).combine(ref_fasta7).combine(ref_fai7)

    output:
    file("${sample_ID}.vcf")
    file("${sample_ID}.norm.vcf")

    script:
    """
    samtools index "${sample_bam}"

    "${params.lofreq_bin}" call-parallel \
    --call-indels \
    --pp-threads \${NSLOTS:-1} \
    --ref "${ref_fasta}" \
    --bed "${targets_bed_file}" \
    --out "${sample_ID}.vcf" \
    "${sample_bam}"

    bgzip -c "${sample_ID}.vcf" > "${sample_ID}.vcf.bgz"

    bcftools index "${sample_ID}.vcf.bgz"

    bcftools norm \
    --multiallelics \
    -both \
    --output-type v \
    "${sample_ID}.vcf.bgz" | \
    bcftools norm \
    --fasta-ref "${ref_fasta}" \
    --output-type v - | \
    bcftools view \
    --exclude 'DP<5' \
    --output-type v >  "${sample_ID}.norm.vcf"

    rm -f "${sample_bam}.bai"
    """
}

process gatk_hc {
    tag { "${sample_ID}" }
    publishDir "${params.wes_output_dir}/vcf_hc", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'
    module 'samtools/1.3'

    input:
    set val(sample_ID), file(sample_bam), file(targets_bed_file), file(ref_fasta), file(ref_fai), file(ref_dict) from samples_dd_ra_rc_bam3.combine(targets_bed7).combine(ref_fasta8).combine(ref_fai8).combine(ref_dict8)

    output:
    file("${sample_ID}.vcf")
    file("${sample_ID}.norm.vcf")

    script:
    """
    samtools index "${sample_bam}"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T HaplotypeCaller \
    -dt NONE \
    --logging_level ERROR \
    -nct \${NSLOTS:-1} \
    --max_alternate_alleles 3 \
    --standard_min_confidence_threshold_for_calling 50 \
    --reference_sequence "${ref_fasta}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${sample_bam}" \
    --out "${sample_ID}.vcf"

    cat "${sample_ID}.vcf" | \
    bcftools norm \
    --multiallelics \
    -both \
    --output-type v - | \
    bcftools norm \
    --fasta-ref "${ref_fasta}" \
    --output-type v - | \
    bcftools view \
    --exclude 'DP<5' \
    --output-type v > "${sample_ID}.norm.vcf"

    rm -f "${sample_bam}.bai"
    """
}
