#!/bin/bash

# submit the Nextflow pipeline as a job to run on the SGE HPC
# NOTE: requires that all nodes can qsub jobs

qsub_logdir="$(readlink -f .)/logs"
mkdir -p "${qsub_logdir}"

job_name="nextflow-exome"

qsub -wd $PWD -o :${qsub_logdir}/ -e :${qsub_logdir}/ -j y -N "$job_name" <<E0F
    set -x
    make exome
E0F