
params.fastq_raw_sheet = "samples.fastq-raw.csv"
params.wes_output_dir = "wes_output"

Channel.fromPath( file(params.fastq_raw_sheet) )
        .splitCsv()
        .map { row ->
            def sample_ID = row[0]
            def read1 = file(row[1])
            def read2 = file(row[2])

            return [ sample_ID, read1, read2 ]
        }
        .into { sample_fastq_raw; sample_fastq_raw2 }

sample_fastq_raw2.println()

process trimmomatic {
    tag { "${sample_ID}-${read1}-${read2}" }
    publishDir "${params.wes_output_dir}/Trimmomatic", mode: 'copy', overwrite: true
    clusterOptions '-pe threaded 1-8'

    input:
    set val(sample_ID), file(read1), file(read2) from sample_fastq_raw

    output:
    file "${sample_ID}-${read1}-${read2}_output.txt"

    script:
    """
    echo "${sample_ID} ${read1} ${read2} \${NSLOTS:-1}" > "${sample_ID}-${read1}-${read2}_output.txt"
    java -Xms16G -Xmx16G -jar ${params.trimmomatic_jar} PE -threads \${NSLOTS:-1} \
    "${read1}" "${read2}" \
    "${sample_ID}_R1.trim.fastq.gz" "${sample_ID}_R1.unpaired.fastq.gz" \
    "${sample_ID}_R2.trim.fastq.gz" "${sample_ID}_R2.unpaired.fastq.gz" \
    ILLUMINACLIP:/ifs/home/id460/ref/contaminants/trimmomatic.fa:2:30:10:1:true TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:35
    """
// java -Xms16G -Xmx16G -jar /local/apps/trimmomatic/0.33/trimmomatic-0.33.jar PE \
// -threads 16 \
// /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/FASTQ-CLEAN/HapMap-B17-1267_R1.fastq.gz \
// /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/FASTQ-CLEAN/HapMap-B17-1267_R2.fastq.gz \
//
// /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/FASTQ-TRIMMED/HapMap-B17-1267_R1.trim.fastq.gz \
// /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/FASTQ-TRIMMED/HapMap-B17-1267_R1.unpaired.fastq.gz \
//
// /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/FASTQ-TRIMMED/HapMap-B17-1267_R2.trim.fastq.gz \
// /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/FASTQ-TRIMMED/HapMap-B17-1267_R2.unpaired.fastq.gz \
// ILLUMINACLIP:/ifs/home/id460/ref/contaminants/trimmomatic.fa:2:30:10:1:true TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:35 \
// 2> /ifs/data/molecpathlab/NGS580_WES-development/sns-demo/logs-trimmomatic/HapMap-B17-1267.txt

}
