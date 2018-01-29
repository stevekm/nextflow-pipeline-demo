
// process delly2 {
//     tag { sample_ID }
//     publishDir "${params.output_dir}/SNV-Delly2", mode: 'move', overwrite: true,
//         saveAs: {filename ->
//             if (filename == 'deletions.vcf') "${sample_ID}_deletions.vcf"
//             else if (filename == 'duplications.vcf') "${sample_ID}_duplications.vcf"
//             else if (filename == 'inversions.vcf') "${sample_ID}_inversions.vcf"
//             else if (filename == 'translocations.vcf') "${sample_ID}_translocations.vcf"
//             else if (filename == 'insertions.vcf') "${sample_ID}_insertions.vcf"
//             else null
//         }
//
//     input:
//     set val(sample_ID), file(bam_gatk_ra_rc), file(bai_gatk_ra_rc) from sample_bam_delly2
//
//     output:
//     file "deletions.vcf"
//     file "duplications.vcf"
//     file "inversions.vcf"
//     file "translocations.vcf"
//     file "insertions.vcf"
//
//     script:
//     """
//     $params.delly2_bin call -t DEL -g ${params.hg19_fa} -o deletions.bcf "${bam_gatk_ra_rc}"
//     $params.delly2_bcftools_bin view deletions.bcf > deletions.vcf
//
//     $params.delly2_bin call -t DUP -g ${params.hg19_fa} -o duplications.bcf "${bam_gatk_ra_rc}"
//     $params.delly2_bcftools_bin view duplications.bcf > duplications.vcf
//
//     $params.delly2_bin call -t INV -g ${params.hg19_fa} -o inversions.bcf "${bam_gatk_ra_rc}"
//     $params.delly2_bcftools_bin view inversions.bcf > inversions.vcf
//
//     $params.delly2_bin call -t BND -g ${params.hg19_fa} -o translocations.bcf "${bam_gatk_ra_rc}"
//     $params.delly2_bcftools_bin view translocations.bcf > translocations.vcf
//
//     $params.delly2_bin call -t INS -g ${params.hg19_fa} -o insertions.bcf "${bam_gatk_ra_rc}"
//     $params.delly2_bcftools_bin view insertions.bcf > insertions.vcf
//     """
// }
