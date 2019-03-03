# snakemake_multiqc 

This is a simple, reproducible example runs fastqc, followed by multiqc in a Snakemake workflow using test fastq files from [ngs-test-data](https://github.com/snakemake-workflows/ngs-test-data).

To run the example run the following from a working directory, and make sure to have Snakemake installed on your system either locally or in a snakemake specific conda environment:

```
git clone --recursive https://github.com/jakevc/snakemake_multiqc.git
cd snakemake_multiqc
snakemake --use-conda
```

The fastqc and multiqc results will be generated in a `data` subdirectory. 

