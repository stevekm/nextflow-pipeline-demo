# Bioinformatics Analysis Pipeline

Analysis pipeline for exome sequencing with tumor-normal paired samples.

__NOTE__: This repo is deprecated. 

- For Nextflow demos, see this repo: https://github.com/stevekm/nextflow-demos

- for Nextflow pipeline developed from this repo, see here: https://github.com/NYU-Molecular-Pathology/NGS580-nf

# Usage

Clone this directory

```bash
git clone https://github.com/stevekm/nextflow-pipeline-demo.git
cd nextflow-pipeline-demo
```

Set up samplesheets

- for the included example data

```bash
make setup
```

- OR for another directory

```bash
./generate-samplesheets.py path/to/fastq_dir
```

Run exome pipeline (recommended to run in `screen`)

```bash
make exome
```
# Included files

- `example-data`: a directory of sample fastq.gz files

- `targets.bed`, `probes.bed`: exome sequencing probes and targets

- `samples.tumor.normal.csv`: tumor-normal sample pairs sheet
