# Bioinformatics Analysis Pipeline

Analysis pipeline for exome sequencing with tumor-normal paired samples.

__NOTE__: Under development

# Usage

Clone this directory

```bash
git clone https://github.com/stevekm/nextflow-pipeline-demo.git
cd nextflow-pipeline-demo
```

Install NextFlow

```bash
make install
```

Set up samplesheets

- for the included example data

```bash
make setup
```

- for another directory

```bash
./generate-samplesheets.py path/to/fastq_dir
```

Run exome pipeline (recommended to run in `screen`)

```bash
make ex
```
