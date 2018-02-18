// Exome sequencing analysis workflow

// pipeline settings; overriden by nextflow.config and CLI args
params.output_dir = "output-exome"
params.bam_dd_ra_rc_gatk_dir = "bam_dd_ra_rc_gatk"



//
// DATA INPUT CHANNELS
//
// targets .bed file
Channel.fromPath( file(params.targets_bed) ).set{ targets_bed }

// reference files
Channel.fromPath( file(params.targets_bed) ).into { targets_bed; targets_bed2; targets_bed3 }
Channel.fromPath( file(params.ref_fa) ).into { ref_fasta; ref_fasta2 }
Channel.fromPath( file(params.ref_fai) ).into { ref_fai; ref_fai2 }
Channel.fromPath( file(params.ref_dict) ).into { ref_dict; ref_dict2 }
Channel.fromPath( file(params.ref_chrom_sizes) ).set{ ref_chrom_sizes }
Channel.fromPath( file(params.trimmomatic_contaminant_fa) ).set{ trimmomatic_contaminant_fa }
Channel.fromPath( file(params.ref_fa_bwa_dir) ).set{ ref_fa_bwa_dir }
Channel.fromPath( file(params.gatk_1000G_phase1_indels_hg19_vcf) ).set{ gatk_1000G_phase1_indels_vcf }
Channel.fromPath( file(params.mills_and_1000G_gold_standard_indels_hg19_vcf) ).set{ mills_and_1000G_gold_standard_indels_vcf }
Channel.fromPath( file(params.dbsnp_ref_vcf) ).set{ dbsnp_ref_vcf }
Channel.fromPath( file(params.cosmic_ref_vcf) ).set{ cosmic_ref_vcf }


// read samples from analysis samplesheet
Channel.fromPath( file(params.samples_analysis_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_ID = row['Sample']
            def reads1 = row['R1'].tokenize( ',' ).collect { file(it) } // comma-sep string into list of files
            def reads2 = row['R2'].tokenize( ',' ).collect { file(it) }
            return [ sample_ID, reads1, reads2 ]
        }
        .into { samples_R1_R2; samples_R1_R2_2 }









//
// PIPELINE TASKS
//

// PREPROCESSING

process fastq_merge {
    // merge the R1 and R2 fastq files into a single fastq each
    tag { "${sample_ID}" }
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    publishDir "${params.output_dir}/fastq-merge", mode: 'copy', overwrite: true

    input:
    set val(sample_ID), file(fastq_r1: "*"), file(fastq_r2: "*") from samples_R1_R2

    output:
    set val(sample_ID), file("${sample_ID}_R1.fastq.gz"), file("${sample_ID}_R2.fastq.gz") into samples_fastq_merged

    script:
    """
    cat ${fastq_r1} > "${sample_ID}_R1.fastq.gz"
    cat ${fastq_r2} > "${sample_ID}_R2.fastq.gz"
    """
}

process trimmomatic {
    // Illumina read trimming
    // http://www.usadellab.org/cms/?page=trimmomatic
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/fastq-trim", mode: 'copy', overwrite: true
    clusterOptions '-pe threaded 1-8'
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"

    input:
    set val(sample_ID), file(read1), file(read2), file(trimmomatic_contaminant_fa) from samples_fastq_merged.combine(trimmomatic_contaminant_fa)

    output:
    set val(sample_ID), file("${sample_ID}_R1.trim.fastq.gz"), file("${sample_ID}_R2.trim.fastq.gz") into samples_fastq_trimmed

    script:
    """
    java -Xms16G -Xmx16G -jar ${params.trimmomatic_jar} PE -threads \${NSLOTS:-1} \
    "${read1}" "${read2}" \
    "${sample_ID}_R1.trim.fastq.gz" "${sample_ID}_R1.unpaired.fastq.gz" \
    "${sample_ID}_R2.trim.fastq.gz" "${sample_ID}_R2.unpaired.fastq.gz" \
    ILLUMINACLIP:${trimmomatic_contaminant_fa}:2:30:10:1:true TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:35
    """
}

process bwa_mem {
    // first pass alignment with BWA
    tag { "${sample_ID}" }
    clusterOptions '-pe threaded 1-8'
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    module 'bwa/0.7.17'

    input:
    set val(sample_ID), file(fastq_R1_trim), file(fastq_R2_trim), file(ref_fa_bwa_dir) from samples_fastq_trimmed.combine(ref_fa_bwa_dir)

    output:
    set val(sample_ID), file("${sample_ID}.sam") into samples_bwa_sam

    script:
    """
    bwa mem -M -v 1 -t \${NSLOTS:-1} -R '@RG\\tID:${sample_ID}\\tSM:${sample_ID}\\tLB:${sample_ID}\\tPL:ILLUMINA' "${ref_fa_bwa_dir}/genome.fa" "${fastq_R1_trim}" "${fastq_R2_trim}" -o "${sample_ID}.sam"
    """
}

process sambamba_view_sort {
    tag { "${sample_ID}" }
    clusterOptions '-pe threaded 1-8'
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"

    input:
    set val(sample_ID), file(sample_sam) from samples_bwa_sam

    output:
    set val(sample_ID), file("${sample_ID}.bam") into samples_bam, samples_bam2

    script:
    """
    "${params.sambamba_bin}" view --sam-input --nthreads=\${NSLOTS:-1} --filter='mapping_quality>=10' --format=bam --compression-level=0 "${sample_sam}" | \
    "${params.sambamba_bin}" sort --nthreads=\${NSLOTS:-1} --memory-limit=16GB --out="${sample_ID}.bam" /dev/stdin
    """
}

process sambamba_flagstat {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/sambamba-flagstat", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"

    input:
    set val(sample_ID), file(sample_bam) from samples_bam

    output:
    file "${sample_ID}.flagstat.txt"

    script:
    """
    "${params.sambamba_bin}" flagstat "${sample_bam}" > "${sample_ID}.flagstat.txt"
    """
}

process sambamba_dedup {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/bam-bwa-dd", mode: 'copy', overwrite: true
    clusterOptions '-pe threaded 1-8'
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    module 'samtools/1.3'

    input:
    set val(sample_ID), file(sample_bam) from samples_bam2

    output:
    set val(sample_ID), file("${sample_ID}.dd.bam") into samples_dd_bam, samples_dd_bam2, samples_dd_bam3, samples_dd_bam4, samples_dd_bam5, samples_dd_bam6, samples_dd_bam7
    file("${sample_ID}.dd.bam.bai")

    script:
    """
    "${params.sambamba_bin}" markdup --remove-duplicates --nthreads \${NSLOTS:-1} --hash-table-size 525000 --overflow-list-size 525000 "${sample_bam}" "${sample_ID}.dd.bam"
    samtools view "${sample_ID}.dd.bam"
    """
}

process sambamba_dedup_flagstat {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/bam-bwa-dd", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"

    input:
    set val(sample_ID), file(sample_bam) from samples_dd_bam2

    output:
    file "${sample_ID}.dd.flagstat.txt"

    script:
    """
    "${params.sambamba_bin}" flagstat "${sample_bam}" > "${sample_ID}.dd.flagstat.txt"
    """

}






// setup downstream Channels
samples_dd_bam.combine(ref_fasta)
            .combine(ref_fai)
            .combine(ref_dict)
            .tap { samples_dd_bam_ref }
            .combine(targets_bed)
            .tap { samples_dd_bam_ref2;
                    samples_dd_bam_ref3;
                    samples_dd_bam_ref4
                }
            .combine(gatk_1000G_phase1_indels_vcf)
            .combine(mills_and_1000G_gold_standard_indels_vcf)
            .combine(dbsnp_ref_vcf)
            .set { samples_dd_bam_ref_gatk }






// GATK RECALIBRATION AND VARIANT CALLING

process qc_target_reads_gatk_genome {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/qc-target-reads", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'

    input:
    set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict) from samples_dd_bam_ref

    output:
    file "${sample_ID}.genome.sample_statistics"
    file "${sample_ID}.genome.sample_summary"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage \
    -dt NONE \
    -rf BadCigar \
    -nt \${NSLOTS:-1} \
    --logging_level ERROR \
    --omitIntervalStatistics \
    --omitLocusTable \
    --omitDepthOutputAtEachBase \
    -ct 10 -ct 100 -mbq 20 -mmq 20 \
    --reference_sequence "${ref_fasta}" \
    --input_file "${sample_bam}" \
    --outputFormat csv \
    --out "${sample_ID}.genome"
    """
}


process qc_target_reads_gatk_pad500 {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/qc-target-reads", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'

    input:
    set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_bam_ref2

    output:
    file "${sample_ID}.pad500.sample_statistics"
    file "${sample_ID}.pad500.sample_summary"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage \
    -dt NONE \
    -rf BadCigar \
    -nt \${NSLOTS:-1} \
    --logging_level ERROR \
    --omitIntervalStatistics \
    --omitLocusTable \
    --omitDepthOutputAtEachBase \
    -ct 10 -ct 100 -mbq 20 -mmq 20 \
    --reference_sequence "${ref_fasta}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 500 \
    --input_file "${sample_bam}" \
    --outputFormat csv \
    --out "${sample_ID}.pad500"
    """
}

process qc_target_reads_gatk_pad100 {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/qc-target-reads", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'

    input:
    set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_bam_ref3

    output:
    file "${sample_ID}.pad100.sample_statistics"
    file "${sample_ID}.pad100.sample_summary"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage \
    -dt NONE \
    -rf BadCigar \
    -nt \${NSLOTS:-1} \
    --logging_level ERROR \
    --omitIntervalStatistics \
    --omitLocusTable \
    --omitDepthOutputAtEachBase \
    -ct 10 -ct 100 -mbq 20 -mmq 20 \
    --reference_sequence "${ref_fasta}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 100 \
    --input_file "${sample_bam}" \
    --outputFormat csv \
    --out "${sample_ID}.pad100"
    """
}

process qc_target_reads_gatk_bed {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/qc-target-reads", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'

    input:
    set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_bam_ref4

    output:
    file "${sample_ID}.bed.sample_statistics"
    file "${sample_ID}.bed.sample_summary"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage \
    -dt NONE \
    -rf BadCigar \
    -nt \${NSLOTS:-1} \
    --logging_level ERROR \
    --omitIntervalStatistics \
    --omitLocusTable \
    --omitDepthOutputAtEachBase \
    -ct 10 -ct 100 -mbq 20 -mmq 20 \
    --reference_sequence "${ref_fasta}" \
    --intervals "${targets_bed_file}" \
    --input_file "${sample_bam}" \
    --outputFormat csv \
    --out "${sample_ID}.bed"
    """
}


// MAIN REALIGNMENT AND RECALIBRATION STEP
process bam_ra_rc_gatk {
    // re-alignment and recalibration with GATK
    // https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
    // https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/${params.bam_dd_ra_rc_gatk_dir}", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'
    module 'samtools/1.3'


    input:
    set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file), file(gatk_1000G_phase1_indels_vcf), file(mills_and_1000G_gold_standard_indels_vcf), file(dbsnp_ref_vcf) from samples_dd_bam_ref_gatk

    output:
    set val(sample_ID), file("${sample_ID}.dd.ra.rc.bam"), file("${sample_ID}.dd.ra.rc.bam.bai") into samples_dd_ra_rc_bam
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
    -known "${gatk_1000G_phase1_indels_vcf}" \
    -known "${mills_and_1000G_gold_standard_indels_vcf}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${sample_bam}" \
    --out "${sample_ID}.intervals"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T IndelRealigner \
    -dt NONE \
    --logging_level ERROR \
    --reference_sequence "${ref_fasta}" \
    --maxReadsForRealignment 50000 \
    -known "${gatk_1000G_phase1_indels_vcf}" \
    -known "${mills_and_1000G_gold_standard_indels_vcf}" \
    -targetIntervals "${sample_ID}.intervals" \
    --input_file "${sample_bam}" \
    --out "${sample_ID}.dd.ra.bam"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T BaseRecalibrator \
    --logging_level ERROR \
    -nct \${NSLOTS:-1} \
    -rf BadCigar \
    --reference_sequence "${ref_fasta}" \
    -knownSites "${gatk_1000G_phase1_indels_vcf}" \
    -knownSites "${mills_and_1000G_gold_standard_indels_vcf}" \
    -knownSites "${dbsnp_ref_vcf}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${sample_ID}.dd.ra.bam" \
    --out "${sample_ID}.table1.txt"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T BaseRecalibrator \
    --logging_level ERROR \
    -nct \${NSLOTS:-1} \
    -rf BadCigar \
    --reference_sequence "${ref_fasta}" \
    -knownSites "${gatk_1000G_phase1_indels_vcf}" \
    -knownSites "${mills_and_1000G_gold_standard_indels_vcf}" \
    -knownSites "${dbsnp_ref_vcf}" \
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

    samtools index "${sample_ID}.dd.ra.rc.bam"
    """
}







// setup downstream Channels
samples_dd_ra_rc_bam.combine(ref_fasta2)
                    .combine(ref_fai2)
                    .combine(ref_dict2)
                    .combine(targets_bed2)
                    .tap { samples_dd_ra_rc_bam_ref;
                            samples_dd_ra_rc_bam_ref2;
                            samples_dd_ra_rc_bam_ref3;
                            samples_dd_ra_rc_bam_ref4;
                            samples_dd_ra_rc_bam_ref5;
                            samples_dd_ra_rc_bam_ref6;
                            samples_dd_ra_rc_bam_ref7 }












process qc_coverage_gatk {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/qc_coverage_gatk", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'

    input:
    set val(sample_ID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref

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
    publishDir "${params.output_dir}/targets", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    module 'bedtools/2.26.0'

    input:
    set file(targets_bed_file), file(ref_chrom_sizes) from targets_bed3.combine(ref_chrom_sizes)

    output:
    file("targets.pad10.bed") into targets_pad_bed

    script:
    """
    cat "${targets_bed_file}" | LC_ALL=C sort -k1,1 -k2,2n | bedtools slop -g "${ref_chrom_sizes}" -b 10 | bedtools merge -d 5 > targets.pad10.bed
    """
}

process lofreq {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/vcf_lofreq", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'
    module 'samtools/1.3'

    input:
    set val(sample_ID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref2

    output:
    file("${sample_ID}.vcf")
    file("${sample_ID}.norm.vcf")

    script:
    """
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
    """
}

process gatk_hc {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/vcf_hc", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'
    module 'samtools/1.3'

    input:
    set val(sample_ID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref3

    output:
    file("${sample_ID}.vcf")
    set val(sample_ID), file("${sample_ID}.norm.vcf") into sample_vcf_hc

    script:
    """
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
    """
}








//
// DOWNSTREAM TASKS
//
// DELLY2 SNV STEPS
process delly2_deletions {
    tag { sample_ID }
    publishDir "${params.output_dir}/snv-deletions-Delly2", mode: 'copy', overwrite: true

    input:
    set val(sample_ID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref4

    output:
    file "${sample_ID}.deletions.vcf"

    script:
    """
    ${params.delly2_bin} call -t DEL -g ${ref_fasta} -o "${sample_ID}.deletions.bcf" "${sample_bam}"
    ${params.delly2_bcftools_bin} view "${sample_ID}.deletions.bcf" > "${sample_ID}.deletions.vcf"
    """
}
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


// Genomic Signatures
process deconstructSigs_signatures {
    tag { "${sample_ID}" }
    validExitStatus 0,11 // allow '11' failure triggered by no variants
    errorStrategy 'ignore'
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    publishDir "${params.output_dir}/signatures_hc", mode: 'copy', overwrite: true,
        saveAs: {filename ->
            if (filename == 'sample_signatures.Rds') "${sample_ID}_signatures.Rds"
            else if (filename == 'sample_signatures.pdf') "${sample_ID}_signatures.pdf"
            else if (filename == 'sample_signatures_pie.pdf') "${sample_ID}_signatures_pie.pdf"
            else filename
        }
    input:
    set val(sample_ID), file(sample_vcf) from sample_vcf_hc

    output:
    file "${sample_ID}_signatures.Rds"
    file "${sample_ID}_signatures.pdf" into signatures_plots
    file "${sample_ID}_signatures_pie.pdf" into signatures_pie_plots

    script:
    """
    deconstructSigs_make_signatures.R "${sample_ID}" "${sample_vcf}"
    """
}



// // REQUIRES ANNOTATIONS FOR DBSNP FILTERING
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


// REQUIRES PAIRED SAMPLES BAM FILES
// process msisensor {
//     tag { sample_ID }
//     module 'samtools/1.3'
//     clusterOptions '-pe threaded 1-8 -j y -l mem_free=40G'
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
//     file "${sample_tumor_ID}_${sample_normal_ID}.msisensor"
//     file "${sample_tumor_ID}_${sample_normal_ID}.msisensor_dis"
//     file "${sample_tumor_ID}_${sample_normal_ID}.msisensor_germline"
//     file "${sample_tumor_ID}_${sample_normal_ID}.msisensor_somatic"
//
//     script:
//     // $params.subset_msisensor_microsatellite_list_script *.bam -m $params.microsatellites -o microsatellites_subset.txt
//     // $params.subset_bam_bed_script *.bam -b $regions_bed -o targets_subset.bed
//     // rm -f "${sample_normal_bam}.bai" "${sample_tumor_bam}.bai" targets_subset.bed microsatellites_subset.txt
//     """
//     samtools index "$sample_tumor_bam"
//     samtools index "$sample_normal_bam"
//     $params.msisensor_bin msi -d $params.microsatellites -n $sample_normal_bam -t $sample_tumor_bam -e $regions_bed -o "${sample_tumor_ID}_${sample_normal_ID}.msisensor" -l 1 -q 1 -b \${NSLOTS:-1}
//     rm -f "${sample_normal_bam}.bai" "${sample_tumor_bam}.bai"
//     """
// }
