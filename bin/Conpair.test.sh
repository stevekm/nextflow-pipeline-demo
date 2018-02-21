#!/bin/bash

module purge
module load python/2.7.3
module load java/1.8
module load jre/1.8
export PYTHONPATH= && \
. activate

export CONPAIR_DIR=/ifs/data/molecpathlab/NGS580_WES-development/nextflow-pipeline-dev/test/Conpair
export GATK_JAR=/ifs/data/sequence/results/external/NYU/snuderllab/bin/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar
export PYTHONPATH="/ifs/data/molecpathlab/NGS580_WES-development/nextflow-pipeline-dev/test/Conpair/modules:${PYTHONPATH}"

results_dir="/ifs/data/molecpathlab/NGS580_WES/170918_NB501073_0025_AHH35JBGX3/results_2017-09-21_15-24-16"


set -x
TUMOR_bam="${results_dir}/BAM-GATK-RA-RC/....dd.ra.rc.bam"
NORMAL_bam="${results_dir}/BAM-GATK-RA-RC/....dd.ra.rc.bam"

"${CONPAIR_DIR}"/scripts/run_gatk_pileup_for_sample.py -B "${TUMOR_bam}" -O TUMOR_pileup
# "${CONPAIR_DIR}"/scripts/run_gatk_pileup_for_sample.py -B "${NORMAL_bam}" -O NORMAL_pileup
#
# "${CONPAIR_DIR}"/scripts/verify_concordance.py -T TUMOR_pileup -N NORMAL_pileup
#
# "${CONPAIR_DIR}"/scripts/estimate_tumor_normal_contamination.py -T TUMOR_pileup -N NORMAL_pileup
