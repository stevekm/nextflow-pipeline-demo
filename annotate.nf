
// Get list of sampleIDs from the raw fastq sheet; map unique sample ID to .vcf file
Channel.fromPath( file(params.samples_fastq_raw_sheet) )
                .splitCsv(sep: ',')
                .map { row ->
                    def sample_ID = row[0]
                    return sample_ID
                }
                .unique()
                .map { sample_ID ->
                    def sample_vcf = file("${params.variants_gatk_hc_dir}/${sample_ID}.original.vcf")
                    return [ sample_ID, sample_vcf ]
                }
                .set{ sample_vcfs }

process annotate_ANNOVAR {
    tag { sample_ID }
    publishDir "${params.output_dir}/Annotations", mode: 'copy', overwrite: true
    clusterOptions '-pe threaded 1-4'

    input:
    set val(sample_ID), file(sample_vcf) from sample_vcfs

    output:
    file "${sample_ID}.${params.build_version}_multianno.txt"
    file "${sample_ID}.sample.${params.build_version}_multianno.txt" into sample_vcfs_summary

    script:
    """
    "${params.ANNOVAR_DIR}/convert2annovar.pl" --format vcf4old --includeinfo "${sample_vcf}" --outfile "${sample_ID}.avinput"

    [ ! "$(wc -l "${sample_ID}.avinput" )" -eq "$(grep -v '^#' "${sample_vcf}" | wc -l)" ] && echo "ERROR: number of lines does not match" && exit 1 || :

    "${params.ANNOVAR_DIR}/table_annovar.pl" "${sample_ID}.avinput" "${params.ANNOVAR_DB_DIR}" --buildver "${params.build_version}" --remove --protocol "refGene,1000g2015aug_all,clinvar_20170905,intervar_20170202,dbnsfp33a,esp6500siv2_all,kaviar_20150923,gnomad_exome,gnomad_genome,avsnp150,fathmm,eigen" --operation "g,f,f,f,f,f,f,f,f,f,f,f" --nastring . --outfile "${sample_ID}" --thread \${NSLOTS:-1}

    paste_col.py -i "${sample_ID}.${params.build_version}_multianno.txt" -o "${sample_ID}.sample.${params.build_version}_multianno.txt" --header "Sample" -v "${sample_ID}" -d "\t"
    """
}

process annotation_summary {
    publishDir "${params.output_dir}/Annotations-Summary", mode: 'copy', overwrite: true

    input:
    file "*" from sample_vcfs_summary.toList()

    script:
    """
    num_lines="$(find . -type f -exec wc -l {} \; | add_col.py)"

    concat_tables.py * > annotation.summary.tsv

    [ ! "$(cat annotation.summary.tsv | wc -l)" -eq "$num_lines" ] && echo "ERROR: number of lines does not match" && exit 1 || :

    """
}
