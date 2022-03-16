import pandas as pd
import os
metadata = pd.read_table("Config/metadata.tsv").set_index("Sample", drop=False)
IDS = metadata["Sample"].tolist()
print(IDS)
os.system("export PATH=$PATH:/mnt/TeacherFiles/Sander_Ferry/Software/signalp-5.0b/bin")


report: "Report/RunSummary.rst"


rule all:
    input:
        expand("QC/raw_reads/{id}", id=IDS),
        expand("QC/filtered_reads/{id}", id=IDS),
        expand("QC/quast/{id}/basic_stats/{id}.draft_GC_content_plot.pdf", id=IDS),
        expand("Funannotate/{id}/logfiles/funannotate-annotate.log", id=IDS),
        "QC/MultiQC/multiqc_report.html"    
    shell:
        "echo {input}"

rule NanoPlot_raw:
    output:
        dir = directory("QC/raw_reads/{sample}"),
        report = report("QC/raw_reads/{sample}/LengthvsQualityScatterPlot_dot.png", caption="Report/NPraw.rst", \
            category="Read quality", subcategory="Raw Reads")
    input:
        "input/{sample}.fastq"
    conda:
        "envs/nanoplot.yaml"
    benchmark:
        "Benchmarks/Nanoplot/{sample}_raw.txt"
    shell:
        "NanoPlot -o {output.dir} --fastq {input}"

rule NanoPlot_trim:
    output:
        dir = directory("QC/filtered_reads/{sample}"),
        report = report("QC/filtered_reads/{sample}/LengthvsQualityScatterPlot_dot.png", caption="Report/NPtrim.rst", \
            category="Read quality", subcategory="Filtered Reads")
    input:
        "filtered_reads/{sample}.fastq"
    conda:
        "envs/nanoplot.yaml"
    benchmark:
        "Benchmarks/Nanoplot/{sample}_filtered.txt"
    shell:
        "NanoPlot -o {output.dir} --fastq {input}"

rule NanoFilt:
    output:
        "filtered_reads/{sample}.fastq"
    input:
        "input/{sample}.fastq"
    conda:
        "envs/preprocessing.yaml"
    benchmark:
        "Benchmarks/NanoFilt/{sample}.txt"
    shell:
        "NanoFilt -l 500 -q 10 < {input} > {output}"

rule Flye:
    output:
        "AssemblyWD/{sample}/assembly.fasta"
    input:
        "filtered_reads/{sample}.fastq"
    conda:
        "envs/preprocessing.yaml"
    threads: 2
    benchmark:
        "Benchmarks/Flye/{sample}.txt"
    shell:
        "flye -g 12.5m -t {threads} -o AssemblyWD/{wildcards.sample} \
        --nano-raw {input} --scaffold"

rule collect_assemblies:
    output:
        "Assemblies/{sample}.draft.fasta"
    input:
        "AssemblyWD/{sample}/assembly.fasta"
    shell:
        "cp {input} {output}"


rule BUSCO:
    output:
        dir = directory("QC/BUSCO/{sample}"),
        sum = "QC/BUSCO/{sample}/short_summary.specific.saccharomycetales_odb10.{sample}.txt"
    input: 
        "Assemblies/{sample}.draft.fasta"
    conda:
        "envs/BUSCO.yaml"
    benchmark:
        "Benchmarks/BUSCO/{sample}.txt"
    shell:
        "busco -i {input} --auto-lineage-euk --out_path QC/BUSCO/ -o {output.dir} -m geno "

rule MultiQC:
    output:
        rep = report("QC/MultiQC/multiqc_report.html", caption = "Report/MQC.rst", category = "Assembly stats", )
    input:
        expand("QC/BUSCO/{sample}/", sample = IDS),
        expand("QC/quast/{sample}/", sample = IDS),
    conda:
        "envs/MultiQC.yaml"
    benchmark:
        "Benchmarks/MultiQC/{sample}.txt"
    shell:
        "multiqc -o {output.dir} QC/"

rule quast:
    output:
        dir = directory("QC/quast/{sample}"),
        report = report("QC/quast/{sample}/basic_stats/{sample}.draft_GC_content_plot.pdf",\
        caption = "Report/quast_GC.rst", category = "Assembly stats", \
        subcategory = "Quast"),
    input:
        "Assemblies/{sample}.draft.fasta"
    conda:
        "envs/preprocessing.yaml"
    benchmark:
        "Benchmarks/Quast/{sample}.txt"
    shell:
        "quast -o {output.dir} {input}"

rule funannotate_clean:
    output:
        "Funannotate/{sample}.clean.fasta"
    input:
        "Assemblies/{sample}.draft.fasta"
    conda:
        "envs/funannotate.yaml"
    benchmark:
        "Benchmarks/Funannotate/{sample}_clean.txt"
    shell:
        "funannotate clean -i {input} -o {output}"

rule funannotate_sort:
    output:
        "Funannotate/{sample}.sort.fasta"
    input:
        "Funannotate/{sample}.clean.fasta"
    conda:
        "envs/funannotate.yaml"
    benchmark:
        "Benchmarks/Funannotate/{sample}_sort.txt"
    shell:
        "funannotate sort -i {input} -o {output} -b contig"

rule funannotate_mask:
    output:
        "Funannotate/{sample}.mask.fasta"
    input:
        "Funannotate/{sample}.sort.fasta"
    conda:
        "envs/funannotate.yaml"
    benchmark:
        "Benchmarks/Funannotate/{sample}_mask.txt"
    shell:
        "funannotate mask -i {input} -o {output}"

rule funannotate_predict:
    output:
        "Funannotate/{sample}/logfiles/funannotate-predict.log"
    input:
        "Funannotate/{sample}.mask.fasta"
    conda:
        "envs/funannotate.yaml"
    benchmark:
        "Benchmarks/Funannotate/{sample}_predict.txt"
    params:
        species = lambda wildcards: metadata.loc[wildcards.sample, "Taxon"],
        GM_path = "/mnt/TeacherFiles/Sander_Ferry/Software/gmes_linux_64",
        db = "/mnt/TeacherFiles/Sander_Ferry/Software/funannotate_db"
    threads: 2
    shell:
        "funannotate predict -i {input} -o Funannotate/{wildcards.sample} \
        -s '{params.species}' --isolate {wildcards.sample} --cpus {threads} \
        --GENEMARK_PATH {params.GM_path} -d {params.db}"

rule funannotate_iprscan:
    output:
        "Funannotate/{sample}/annotate_misc/iprscan.xml"
    input:
        "Funannotate/{sample}/logfiles/funannotate-predict.log"
    conda:
        "envs/funannotate.yaml"
    benchmark:
        "Benchmarks/Funannotate/{sample}_IPR.txt"
    params:
        ipr = "/mnt/TeacherFiles/Sander_Ferry/Software/interproscan-5.52-86.0/"
    threads: 2
    shell:
        "funannotate iprscan -i Funannotate/{wildcards.sample} -m local \
        --cpus {threads} --iprscan_path {params.ipr}"


rule funannotate_annotate:
    output:
        "Funannotate/{sample}/logfiles/funannotate-annotate.log"
    input:
        log = "Funannotate/{sample}/logfiles/funannotate-predict.log",
        ipr = "Funannotate/{sample}/annotate_misc/iprscan.xml"
    conda:
        "envs/funannotate.yaml"
    benchmark:
        "Benchmarks/Funannotate/{sample}_annotate.txt"
    params:
        sbt = "Config/NCBI.sbt",
        species = lambda wildcards: metadata.loc[wildcards.sample, "Taxon"],
        db = "/mnt/TeacherFiles/Sander_Ferry/Software/funannotate_db"
    threads: 4
    shell:
        "funannotate annotate -i Funannotate/{wildcards.sample} \
        --sbt {params.sbt} --species '{params.species}' \
        --iprscan {input.ipr} -d {params.db} --isolate {wildcards.sample}"
