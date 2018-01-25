//
// pipeline to run MSI Sensor on a single sample
// https://github.com/ding-lab/msisensor
//

params.bam_dir = "BAM-BWA"
params.samples_pairs_sheet = "samples.pairs.controls.csv"
params.regions_bed = "targets.bed"
params.msisensor_bin = "/ifs/data/molecpathlab/bin/msisensor/msisensor"
params.microsatellites = "/ifs/data/molecpathlab/bin/msisensor/hg19_microsatellites.list"


// params.sampleID = null
// params.bam_normal = null
// params.bam_tumor = null
// check mandatory options
// if (!params.sampleID) {
//     exit 1, "Sample ID genome not specified"
// }
// if (!params.bam_normal) {
//     exit 1, "Sample normal .bam file not specified"
// }
// if (!params.bam_tumor) {
//     exit 1, "Sample tumor .bam file not specified"
// }
// println """\
//          MSI Sensor
//          ===================================
//          sampleID: ${params.sampleID}
//          bam_normal        : ${params.bam_normal}
//          bam_tumor        : ${params.bam_tumor}
//          regions_bed        : ${params.regions_bed}
//          msisensor_bin       : ${params.msisensor_bin}
//          microsatellites       : ${params.microsatellites}
//          """
//          .stripIndent()


/*
* create file objects from the given string parameters
*/
bam_normal_file = file(params.bam_normal)
bai_normal_file = file("${params.bam_normal}.bai")
bam_tumor_file = file(params.bam_tumor)
bai_tumor_file = file("${params.bam_tumor}.bai")
microsatellites_file = file(params.microsatellites)
regions_bed_file = file(params.regions_bed)

process msisensor {
    // run msisensor
    clusterOptions '-cwd -pe threaded 1-32 -j y'
    publishDir "output/${params.sampleID}/MSISensor" , mode: 'copy', overwrite: true,
        saveAs: {filename ->
            if (filename == 'msisensor') "${params.sampleID}.msisensor"
            else if (filename == 'msisensor_dis') "${params.sampleID}.msisensor_dis"
            else if (filename == 'msisensor_germline') "${params.sampleID}.msisensor_germline"
            else if (filename == 'msisensor_somatic') "${params.sampleID}.msisensor_somatic"
            else null
        }

    input:
    file microsatellites from microsatellites_file
    file bam_normal from bam_normal_file
    file bai_normal from bai_normal_file
    file bam_tumor from bam_tumor_file
    file bai_tumor from bai_tumor_file
    file regions_bed from regions_bed_file

    output:
    file 'msisensor'
    file 'msisensor_dis'
    file 'msisensor_germline'
    file 'msisensor_somatic'

    script:
    """
    $params.msisensor_bin msi -d $microsatellites -n $bam_normal -t $bam_tumor -e $regions_bed -o msisensor -l 1 -q 1 -b \${NSLOTS:-1}
    """
}
