params.fastq_raw_sheet = "samples.fastq-raw.csv"
params.targets_bed = "targets.bed"
params.wes_output_dir = "output-exomes"

// read samples from fastq raw sheet
Channel.fromPath( file(params.fastq_raw_sheet) )
        .splitCsv()
        .map { row ->
            def sample_ID = row[0]
            return sample_ID
        }
        .set { sample_IDs }


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

Channel.fromPath( file(params.targets_bed) )
        .into { targets_bed; targets_bed2; targets_bed3; targets_bed4 }
//
//
// DEBUGGING
//
//

// print the unique sample IDs
sample_IDs.unique().println()

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

process fastq_merge {
    // merge the R1 and R2 fastq files into a single fastq each
    tag { "${sample_ID}" }
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    publishDir "${params.wes_output_dir}/fastq-merge", mode: 'copy', overwrite: true

    input:
    set val(sample_ID), file(fastq_r1: "*"), file(fastq_r2: "*") from sample_fastq_r1r2

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
    tag { "${sample_ID}-${read1}-${read2}" }
    publishDir "${params.wes_output_dir}/fastq-trim", mode: 'copy', overwrite: true
    clusterOptions '-pe threaded 1-8'
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"

    input:
    set val(sample_ID), file(read1), file(read2) from samples_fastq_merged

    output:
    set val(sample_ID), file("${sample_ID}_R1.trim.fastq.gz"), file("${sample_ID}_R1.unpaired.fastq.gz"), file("${sample_ID}_R2.trim.fastq.gz"), file("${sample_ID}_R2.unpaired.fastq.gz") into samples_fastq_trimmed

    script:
    """
    # echo "${sample_ID} ${read1} ${read2} \${NSLOTS:-1}"
    java -Xms16G -Xmx16G -jar ${params.trimmomatic_jar} PE -threads \${NSLOTS:-1} \
    "${read1}" "${read2}" \
    "${sample_ID}_R1.trim.fastq.gz" "${sample_ID}_R1.unpaired.fastq.gz" \
    "${sample_ID}_R2.trim.fastq.gz" "${sample_ID}_R2.unpaired.fastq.gz" \
    ILLUMINACLIP:${params.trimmomatic_contaminant_fa}:2:30:10:1:true TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:35
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
    set val(sample_ID), file(fastq_R1_trim), file(fastq_R1_trim_unpaired), file(fastq_R2_trim), file(fastq_R2_trim_unpaired) from samples_fastq_trimmed

    output:
    set val(sample_ID), file("${sample_ID}.sam") into samples_bwa_sam

    script:
    """
    bwa mem -M -v 1 -t \${NSLOTS:-1} -R '@RG\\tID:${sample_ID}\\tSM:${sample_ID}\\tLB:${sample_ID}\\tPL:ILLUMINA' "${params.bwa_hg19_ref_fa}" "${fastq_R1_trim}" "${fastq_R2_trim}" -o "${sample_ID}.sam"
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
    publishDir "${params.wes_output_dir}/sambamba-flagstat", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"

    input:
    set val(sample_ID), file(sample_bam) from samples_bam2

    output:
    file "${sample_ID}.flagstat.txt"

    script:
    """
    "${params.sambamba_bin}" flagstat "${sample_bam}" > "${sample_ID}.flagstat.txt"
    """
}

process sambamba_dedup {
    tag { "${sample_ID}" }
    publishDir "${params.wes_output_dir}/bam-bwa", mode: 'copy', overwrite: true
    clusterOptions '-pe threaded 1-8'
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"

    input:
    set val(sample_ID), file(sample_bam) from samples_bam

    output:
    set val(sample_ID), file("${sample_ID}.dd.bam") into samples_dd_bam, samples_dd_bam2, samples_dd_bam3, samples_dd_bam4, samples_dd_bam5, samples_dd_bam6

    script:
    """
    "${params.sambamba_bin}" markdup --remove-duplicates --nthreads \${NSLOTS:-1} --hash-table-size 525000 --overflow-list-size 525000 "${sample_bam}" "${sample_ID}.dd.bam"
    """
}

process sambamba_dedup_flagstat {
    tag { "${sample_ID}" }
    // executor "local"
    publishDir "${params.wes_output_dir}/sambamba-flagstat-dd", mode: 'copy', overwrite: true
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

process qc_target_reads_gatk_genome {
    tag { "${sample_ID}" }
    // executor "local"
    publishDir "${params.wes_output_dir}/qc-target-reads", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'

    input:
    set val(sample_ID), file(sample_bam) from samples_dd_bam

    output:
    file "${sample_ID}.genome.sample_statistics"
    file "${sample_ID}.genome.sample_summary"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage -dt NONE -rf BadCigar -nt \${NSLOTS:-1} --logging_level ERROR --omitIntervalStatistics --omitLocusTable --omitDepthOutputAtEachBase -ct 10 -ct 100 -mbq 20 -mmq 20 --reference_sequence "${params.hg19_fa}" --input_file "${sample_bam}" --outputFormat csv --out "${sample_ID}.genome"
    """
    // java -Xms16G -Xmx16G -jar /ifs/home/id460/software/GenomeAnalysisTK/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T DepthOfCoverage -dt NONE -rf BadCigar -nt 16 --logging_level ERROR --omitIntervalStatistics --omitLocusTable --omitDepthOutputAtEachBase -ct 10 -ct 100 -mbq 20 -mmq 20 --reference_sequence /ifs/data/sequence/Illumina/igor/ref/hg19/genome.fa --input_file /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/BAM-DD/HapMap-B17-1267.dd.bam --outputFormat csv --out /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/QC-target-reads/HapMap-B17-1267.genome

}

process qc_target_reads_gatk_pad500 {
    tag { "${sample_ID}" }
    // executor "local"
    publishDir "${params.wes_output_dir}/qc-target-reads", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'

    input:
    set val(sample_ID), file(sample_bam) from samples_dd_bam3
    file targets_bed_file from targets_bed

    output:
    file "${sample_ID}.pad500.sample_statistics"
    file "${sample_ID}.pad500.sample_summary"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage -dt NONE -rf BadCigar -nt \${NSLOTS:-1} --logging_level ERROR --omitIntervalStatistics --omitLocusTable --omitDepthOutputAtEachBase -ct 10 -ct 100 -mbq 20 -mmq 20 --reference_sequence "${params.hg19_fa}" --intervals "${targets_bed_file}" --interval_padding 500 --input_file "${sample_bam}" --outputFormat csv --out "${sample_ID}.pad500"
    """
    // java -Xms16G -Xmx16G -jar /ifs/home/id460/software/GenomeAnalysisTK/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T DepthOfCoverage -dt NONE -rf BadCigar -nt 16 --logging_level ERROR --omitIntervalStatistics --omitLocusTable --omitDepthOutputAtEachBase -ct 10 -ct 100 -mbq 20 -mmq 20 --reference_sequence /ifs/data/sequence/Illumina/igor/ref/hg19/genome.fa --intervals /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/targets.bed --interval_padding 500 --input_file /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/BAM-DD/HapMap-B17-1267.dd.bam --outputFormat csv --out /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/QC-target-reads/HapMap-B17-1267.pad500
}

process qc_target_reads_gatk_pad100 {
    tag { "${sample_ID}" }
    publishDir "${params.wes_output_dir}/qc-target-reads", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'

    input:
    set val(sample_ID), file(sample_bam) from samples_dd_bam4
    file targets_bed_file from targets_bed2

    output:
    file "${sample_ID}.pad100.sample_statistics"
    file "${sample_ID}.pad100.sample_summary"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage -dt NONE -rf BadCigar -nt \${NSLOTS:-1} --logging_level ERROR --omitIntervalStatistics --omitLocusTable --omitDepthOutputAtEachBase -ct 10 -ct 100 -mbq 20 -mmq 20 --reference_sequence "${params.hg19_fa}" --intervals "${targets_bed_file}" --interval_padding 100 --input_file "${sample_bam}" --outputFormat csv --out "${sample_ID}.pad100"
    """
    // java -Xms16G -Xmx16G -jar /ifs/home/id460/software/GenomeAnalysisTK/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T DepthOfCoverage -dt NONE -rf BadCigar -nt 16 --logging_level ERROR --omitIntervalStatistics --omitLocusTable --omitDepthOutputAtEachBase -ct 10 -ct 100 -mbq 20 -mmq 20 --reference_sequence /ifs/data/sequence/Illumina/igor/ref/hg19/genome.fa --intervals /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/targets.bed --interval_padding 100 --input_file /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/BAM-DD/HapMap-B17-1267.dd.bam --outputFormat csv --out /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/QC-target-reads/HapMap-B17-1267.pad100
}

process qc_target_reads_gatk_bed {
    tag { "${sample_ID}" }
    publishDir "${params.wes_output_dir}/qc-target-reads", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'

    input:
    set val(sample_ID), file(sample_bam) from samples_dd_bam5
    file targets_bed_file from targets_bed3

    output:
    file "${sample_ID}.bed.sample_statistics"
    file "${sample_ID}.bed.sample_summary"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage -dt NONE -rf BadCigar -nt \${NSLOTS:-1} --logging_level ERROR --omitIntervalStatistics --omitLocusTable --omitDepthOutputAtEachBase -ct 10 -ct 100 -mbq 20 -mmq 20 --reference_sequence "${params.hg19_fa}" --intervals "${targets_bed_file}" --input_file "${sample_bam}" --outputFormat csv --out "${sample_ID}.bed"
    """
    // java -Xms16G -Xmx16G -jar /ifs/home/id460/software/GenomeAnalysisTK/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T DepthOfCoverage -dt NONE -rf BadCigar -nt 16 --logging_level ERROR --omitIntervalStatistics --omitLocusTable --omitDepthOutputAtEachBase -ct 10 -ct 100 -mbq 20 -mmq 20 --reference_sequence /ifs/data/sequence/Illumina/igor/ref/hg19/genome.fa --intervals /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/targets.bed --input_file /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/BAM-DD/HapMap-B17-1267.dd.bam --outputFormat csv --out /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/QC-target-reads/HapMap-B17-1267.bed
}

process bam_ra_rc_gatk {
    // re-alignment and recalibration with GATK
    // https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
    // https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php
    tag { "${sample_ID}" }
    publishDir "${params.wes_output_dir}/Alignments", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-8'

    input:
    set val(sample_ID), file(sample_bam) from samples_dd_bam6
    file targets_bed_file from targets_bed4

    output:
    file "${sample_ID}.intervals"
    file "${sample_ID}.dd.ra.rc.bam"
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
    --reference_sequence "${params.hg19_fa}" \
    -known "${params.gatk_1000G_phase1_indels_hg19_vcf}" \
    -known "${params.mills_and_1000G_gold_standard_indels_hg19_vcf}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${sample_bam}" \
    --out "${sample_ID}.intervals"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T IndelRealigner \
    -dt NONE \
    --logging_level ERROR \
    --reference_sequence "${params.hg19_fa}" \
    --maxReadsForRealignment 50000 \
    -known "${params.gatk_1000G_phase1_indels_hg19_vcf}" \
    -known "${params.mills_and_1000G_gold_standard_indels_hg19_vcf}" \
    -targetIntervals "${sample_ID}.intervals" \
    --input_file "${sample_bam}" \
    --out "${sample_ID}.dd.ra.bam"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T BaseRecalibrator \
    --logging_level ERROR \
    -nt \${NSLOTS:-1} \
    -rf BadCigar \
    --reference_sequence "${params.hg19_fa}" \
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
    --reference_sequence "${params.hg19_fa}" \
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
    --reference_sequence "${params.hg19_fa}" \
    -before "${sample_ID}.table1.txt" \
    -after "${sample_ID}.table2.txt" \
    -csv "${sample_ID}.csv" \
    -plots "${sample_ID}.pdf"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T PrintReads \
    --logging_level ERROR \
    -nct \${NSLOTS:-1} \
    -rf BadCigar \
    --reference_sequence "${params.hg19_fa}" \
    -BQSR "${sample_ID}.table1.txt" \
    --input_file "${sample_ID}.dd.ra.bam" \
    --out "${sample_ID}.dd.ra.rc.bam"
    """

// java -Xms16G -Xmx16G -jar /ifs/home/id460/software/GenomeAnalysisTK/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T RealignerTargetCreator -dt NONE --logging_level ERROR -nt 16 --reference_sequence /ifs/data/sequence/Illumina/igor/ref/hg19/genome.fa -known /ifs/home/id460/ref/hg19/gatk-bundle/1000G_phase1.indels.hg19.vcf -known /ifs/home/id460/ref/hg19/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf --intervals /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/targets.bed --interval_padding 10 --input_file /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/BAM-DD/HapMap-B17-1267.dd.bam --out /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/logs-bam-ra-rc-gatk/HapMap-B17-1267.intervals

// java -Xms16G -Xmx16G -jar /ifs/home/id460/software/GenomeAnalysisTK/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T IndelRealigner -dt NONE --logging_level ERROR --reference_sequence /ifs/data/sequence/Illumina/igor/ref/hg19/genome.fa --maxReadsForRealignment 50000 -known /ifs/home/id460/ref/hg19/gatk-bundle/1000G_phase1.indels.hg19.vcf -known /ifs/home/id460/ref/hg19/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf -targetIntervals /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/logs-bam-ra-rc-gatk/HapMap-B17-1267.intervals --input_file /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/BAM-DD/HapMap-B17-1267.dd.bam --out /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/BAM-GATK-RA/HapMap-B17-1267.dd.ra.bam

// java -Xms16G -Xmx16G -jar /ifs/home/id460/software/GenomeAnalysisTK/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T BaseRecalibrator --logging_level ERROR -nct 16 -rf BadCigar --reference_sequence /ifs/data/sequence/Illumina/igor/ref/hg19/genome.fa -knownSites /ifs/home/id460/ref/hg19/gatk-bundle/1000G_phase1.indels.hg19.vcf -knownSites /ifs/home/id460/ref/hg19/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf -knownSites /ifs/home/id460/ref/hg19/gatk-bundle/dbsnp_138.hg19.vcf --intervals /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/targets.bed --interval_padding 10 --input_file /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/BAM-GATK-RA/HapMap-B17-1267.dd.ra.bam --out /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/logs-bam-ra-rc-gatk/HapMap-B17-1267.table1.txt

// java -Xms16G -Xmx16G -jar /ifs/home/id460/software/GenomeAnalysisTK/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T BaseRecalibrator --logging_level ERROR -nct 16 -rf BadCigar --reference_sequence /ifs/data/sequence/Illumina/igor/ref/hg19/genome.fa -knownSites /ifs/home/id460/ref/hg19/gatk-bundle/1000G_phase1.indels.hg19.vcf -knownSites /ifs/home/id460/ref/hg19/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf -knownSites /ifs/home/id460/ref/hg19/gatk-bundle/dbsnp_138.hg19.vcf --intervals /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/targets.bed --interval_padding 10 --input_file /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/BAM-GATK-RA/HapMap-B17-1267.dd.ra.bam -BQSR /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/logs-bam-ra-rc-gatk/HapMap-B17-1267.table1.txt --out /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/logs-bam-ra-rc-gatk/HapMap-B17-1267.table2.txt

// java -Xms16G -Xmx16G -jar /ifs/home/id460/software/GenomeAnalysisTK/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T AnalyzeCovariates --logging_level ERROR --reference_sequence /ifs/data/sequence/Illumina/igor/ref/hg19/genome.fa -before /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/logs-bam-ra-rc-gatk/HapMap-B17-1267.table1.txt -after /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/logs-bam-ra-rc-gatk/HapMap-B17-1267.table2.txt -csv /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/logs-bam-ra-rc-gatk/HapMap-B17-1267.csv -plots /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/logs-bam-ra-rc-gatk/HapMap-B17-1267.pdf

java -Xms16G -Xmx16G -jar /ifs/home/id460/software/GenomeAnalysisTK/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T PrintReads --logging_level ERROR -nct 16 -rf BadCigar --reference_sequence /ifs/data/sequence/Illumina/igor/ref/hg19/genome.fa -BQSR /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/logs-bam-ra-rc-gatk/HapMap-B17-1267.table1.txt --input_file /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/BAM-GATK-RA/HapMap-B17-1267.dd.ra.bam --out /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/BAM-GATK-RA-RC/HapMap-B17-1267.dd.ra.rc.bam
}
