# Bioinformatics Analysis Pipeline

## Usage

Clone this directory

```bash
git clone https://github.com/stevekm/nextflow-pipeline-demo.git
cd nextflow-pipeline-demo
```

Install NextFlow

```bash
make install
```

Make a symlink to your sns analysis output location

```bash
ln -fs /path/to/sns_output sns-dir
```

Run (recommended to run in `screen`)

```bash
./nextflow run pipeline.nf
```
