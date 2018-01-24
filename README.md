Install NextFlow

```bash
curl -fsSL get.nextflow.io | bash
```

Run

```bash
./nextflow run msisensor.nf --sampleID foo --bam_normal bar --bam_tumor baz
```

example

```bash
tumor_bam="/ifs/data/molecpathlab/NGS580_WES-development/msisensor-dev/BAM-BWA/SeraCare-1to1-Positive.bam"
normal_bam="/ifs/data/molecpathlab/NGS580_WES-development/msisensor-dev/BAM-BWA/HapMap-B17-1267.bam"
./nextflow run msisensor.nf --sampleID test1 --bam_normal "$normal_bam" --bam_tumor "$tumor_bam"
```
