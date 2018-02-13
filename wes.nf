
params.fastq_raw_sheet = "samples.fastq-raw.test.csv"
params.wes_output_dir = "wes_output"

// Channel.fromPath( file(params.fastq_raw_sheet) )
//         .splitCsv()
//         .map { row ->
//             def sample_ID = row[0]
//             return sample_ID
//         }
//         .set { sample_IDs }

// sample_IDs.unique().println()

// Channel.fromPath( file(params.fastq_raw_sheet) )
//         .splitCsv()
//         .map { row ->
//             def sample_ID = row[0]
//             def read1 = file(row[1])
//             def read2 = file(row[2])
//
//             return [ sample_ID, read1, read2 ]
//         }
//         .into { sample_fastq_raw; sample_fastq_raw2 }

// sample_fastq_raw2.println()

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

process fastq_pairs_print {
    executor "local"
    echo true
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

process fastq_merge {
    tag { "${sample_ID}" }
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
    tag { "${sample_ID}-${read1}-${read2}" }
    publishDir "${params.wes_output_dir}/fastq-trim", mode: 'copy', overwrite: true
    clusterOptions '-pe threaded 1-8'
    beforeScript 'echo "USER:\${USER:-none}\tJOB_ID:\${JOB_ID:-none}\tJOB_NAME:\${JOB_NAME:-none}\tHOSTNAME:\${HOSTNAME:-none}\t"'

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
    tag { "${sample_ID}" }
    publishDir "${params.wes_output_dir}/bam-bwa", mode: 'copy', overwrite: true
    clusterOptions '-pe threaded 4-16'
    beforeScript 'echo "USER:\${USER:-none}\tJOB_ID:\${JOB_ID:-none}\tJOB_NAME:\${JOB_NAME:-none}\tHOSTNAME:\${HOSTNAME:-none}\t"'
    module 'bwa/0.7.17'

    input:
    set val(sample_ID), file(fastq_R1_trim), file(fastq_R1_trim_unpaired), file(fastq_R2_trim), file(fastq_R2_trim_unpaired) from samples_fastq_trimmed

    output:
    set val(sample_ID), file("${sample_ID}.sam") into samples_bwa_sam

    script:
    """
    bwa mem -M -v 1 -t 16 -R '@RG\\tID:${sample_ID}\\tSM:${sample_ID}\\tLB:${sample_ID}\\tPL:ILLUMINA' "${params.bwa_hg19_ref_fa}" "${fastq_R1_trim}" "${fastq_R2_trim}" -o "${sample_ID}.sam"
    """
    // bwa mem -M -v 1 -t 16 -R '@RG\tID:HapMap-B17-1267\tSM:HapMap-B17-1267\tLB:HapMap-B17-1267\tPL:ILLUMINA' /ifs/data/sequence/Illumina/igor/ref/hg19/BWA/genome.fa /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/FASTQ-TRIMMED/HapMap-B17-1267_R1.trim.fastq.gz /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/FASTQ-TRIMMED/HapMap-B17-1267_R2.trim.fastq.gz | /ifs/home/id460/software/sambamba/sambamba_v0.6.7 view --sam-input --nthreads=16 --filter='mapping_quality>=10' --format=bam --compression-level=0 /dev/stdin | /ifs/home/id460/software/sambamba/sambamba_v0.6.7 sort --nthreads=16 --memory-limit=16GB --out=/ifs/data/molecpathlab/NGS580_WES-development/sns-demo/BAM-BWA/HapMap-B17-1267.bam /dev/stdin
}

process sambamba {
    tag { "${sample_ID}" }
    publishDir "${params.wes_output_dir}/bam-bwa", mode: 'copy', overwrite: true
    clusterOptions '-pe threaded 4-16'
    beforeScript 'echo "USER:\${USER:-none}\tJOB_ID:\${JOB_ID:-none}\tJOB_NAME:\${JOB_NAME:-none}\tHOSTNAME:\${HOSTNAME:-none}\t"'

    input:
    set val(sample_ID), file(sample_sam) from samples_bwa_sam

    output:
    set val(sample_ID), file("${sample_ID}.bam") into samples_bam

    script:
    """
    "${params.sambamba_bin}" view --sam-input --nthreads=\${NSLOTS:-1} --filter='mapping_quality>=10' --format=bam --compression-level=0 "sample_sam" | \
    "${params.sambamba_bin}" sort --nthreads=\${NSLOTS:-1} --memory-limit=16GB --out="${sample_ID}.bam" /dev/stdin
    """
    /ifs/home/id460/software/sambamba/sambamba_v0.6.7 view --sam-input --nthreads=16 --filter='mapping_quality>=10' --format=bam --compression-level=0 /dev/stdin | /ifs/home/id460/software/sambamba/sambamba_v0.6.7 sort --nthreads=16 --memory-limit=16GB --out=/ifs/data/molecpathlab/NGS580_WES-development/sns-demo/BAM-BWA/HapMap-B17-1267.bam /dev/stdin
}
