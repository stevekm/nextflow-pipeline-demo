params.samples_tumor_normal_sheet = "samples.tumor.normal.csv"

Channel.fromPath( file(params.targets_bed) ).set{ targets_bed }
Channel.fromPath( file(params.hg19_fa) ).set{ ref_fasta }
Channel.fromPath( file(params.hg19_fai) ).set{ ref_fai }
Channel.fromPath( file(params.hg19_dict) ).set{ ref_dict }
Channel.fromPath( file(params.dbsnp_ref_vcf) ).set{ dbsnp_ref_vcf }
Channel.fromPath( file(params.cosmic_ref_vcf) ).set{ cosmic_ref_vcf }

// read samples from fastq raw sheet
Channel.fromPath( file(params.fastq_raw_sheet) )
        .splitCsv()
        .map { row ->
            def sample_ID = row[0]
            return sample_ID
        }
        .set { sample_IDs }

// read sample tumor - normal pairs from sheet; only non-NA entries that match the raw fastq sheet samples
Channel.fromPath( file(params.samples_tumor_normal_sheet) )
        .splitCsv()
        .map { row ->
            def sample_ID_tumor = row[0]
            def sample_ID_normal = row[1]
            return [ sample_ID_tumor, sample_ID_normal ]
        }
        .filter{ sample_ID_tumor, sample_ID_normal ->
                        sample_ID_tumor != "NA" && sample_ID_normal != "NA"
        }
        .join(sample_IDs.unique())
        .map { sample_ID_tumor, sample_ID_normal ->
            def comparison_ID = "${sample_ID_tumor}_${sample_ID_normal}"
            def tumor_bam = file("${params.wes_output_dir}/${params.bam_ra_rc_gatk_dir}/${sample_ID_tumor}.dd.ra.rc.bam")
            def normal_bam = file("${params.wes_output_dir}/${params.bam_ra_rc_gatk_dir}/${sample_ID_normal}.dd.ra.rc.bam")
            return [ comparison_ID, sample_ID_tumor, tumor_bam, sample_ID_normal, normal_bam]
        }
        .combine(targets_bed)
        .combine(ref_fasta)
        .combine(ref_fai)
        .combine(ref_dict)
        .combine(dbsnp_ref_vcf)
        .combine(cosmic_ref_vcf)
        .into { samples_pairs; samples_pairs2 }

// get the unique chromosomes in the targets bed file
Channel.fromPath( params.targets_bed )
            .splitCsv(sep: '\t')
            .map{row ->
                row[0]
            }
            .unique()
            .set{ chroms }


process mutect2 {
    tag { "${comparison_ID}:${chrom}" }
    publishDir "${params.wes_output_dir}/vcf_mutect2", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-l mem_free=150G -hard'
    module 'samtools/1.3'
    module 'java/1.8'

    input:
    set val(chrom), val(comparison_ID), val(sample_ID_tumor), file(tumor_bam), val(sample_ID_normal), file(normal_bam), file(targets_bed), file(ref_fasta), file(ref_fai), file(ref_dict), file(dbsnp_ref_vcf), file(cosmic_ref_vcf) from chroms.combine(samples_pairs)

    output:
    file("${comparison_ID}.vcf")

    script:
    """
    subset_bed.py "${chrom}" "${targets_bed}" > "${comparison_ID}.${chrom}.bed"

    samtools index "${tumor_bam}"
    samtools index "${normal_bam}"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T MuTect2 \
    -dt NONE \
    --logging_level WARN \
    --standard_min_confidence_threshold_for_calling 30 \
    --max_alt_alleles_in_normal_count 10 \
    --max_alt_allele_in_normal_fraction 0.05 \
    --max_alt_alleles_in_normal_qscore_sum 40 \
    --reference_sequence "${ref_fasta}" \
    --dbsnp "${dbsnp_ref_vcf}" \
    --cosmic "${cosmic_ref_vcf}" \
    --intervals "${comparison_ID}.${chrom}.bed" \
    --interval_padding 10 \
    --input_file:tumor "${tumor_bam}" \
    --input_file:normal "${normal_bam}" \
    --out "${comparison_ID}.vcf"

    rm -f "${tumor_bam}.bai"
    rm -f "${normal_bam}.bai"
    rm -f "${comparison_ID}.${chrom}.bed"
    """
}
