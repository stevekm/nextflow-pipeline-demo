// Pre-processing steps for exome analysis
params.output_dir = "output-preprocessing"

// targets .bed file
Channel.fromPath( file(params.targets_bed) ).set{ targets_bed }

// reference hg19 fasta file
Channel.fromPath( file(params.targets_bed) ).set{ targets_bed }
Channel.fromPath( file(params.ref_fa) ).set{ ref_fasta }
Channel.fromPath( file(params.ref_fai) ).set{ ref_fai }
Channel.fromPath( file(params.ref_dict) ).set{ ref_dict }
Channel.fromPath( file(params.ref_chrom_sizes) ).set{ ref_chrom_sizes }
Channel.fromPath( file(params.trimmomatic_contaminant_fa) ).set{ trimmomatic_contaminant_fa }
Channel.fromPath( file(params.ref_fa_bwa_dir) ).set{ ref_fa_bwa_dir }


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
            .tap { samples_qc_target_reads_gatk }
            .combine(targets_bed)
            .into { samples_qc_target_reads_gatk_pad500;
                    samples_qc_target_reads_gatk_pad100;
                    samples_qc_target_reads_gatk_bed
                }







process qc_target_reads_gatk_genome {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/qc-target-reads", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'

    input:
    set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict) from samples_qc_target_reads_gatk

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
    set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_qc_target_reads_gatk_pad500

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
    set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_qc_target_reads_gatk_pad100

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
    set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_qc_target_reads_gatk_bed

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
