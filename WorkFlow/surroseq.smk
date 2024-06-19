configfile: "config/config.yaml"


from util import getSamples, getPair
import pandas as pd

###############################################################
#####################功能函数##################################
###############################################################
SAMPLES = getSamples(config["samples"])
PAIRS = getPair(config["samples"])


def getFastqFiles(wildcards):
    df = pd.read_table(config["samples"])
    df = df[df.SampleName == wildcards.sample].reset_index()
    if df.shape[0] != 1:
        raise Exception("Sample name not unique or not found in samples.tsv")
    fqs = [df.loc[0, "FastQ1"], df.loc[0, "FastQ2"]]
    return fqs


def getPairs(wildcards):
    dct = {}
    dct["Case"] = "extractReads/" + PAIRS[wildcards.sample][0]
    dct["Control"] = "extractReads/" + PAIRS[wildcards.sample][1]
    return dct


def getPairs_for_qc(wildcards):
    dct = {}
    dct["Case"] = "fastQC/" + PAIRS[wildcards.sample][0]
    dct["Control"] = "fastQC/" + PAIRS[wildcards.sample][1]
    return dct


def getPairs_for_overlap(wildcards):
    dct = {}
    dct["Case"] = (
        "flash/"
        + PAIRS[wildcards.sample][0]
        + "/"
        + PAIRS[wildcards.sample][0]
        + ".extendedFrags.fastq.gz"
    )
    dct["Control"] = (
        "flash/"
        + PAIRS[wildcards.sample][1]
        + "/"
        + PAIRS[wildcards.sample][1]
        + ".extendedFrags.fastq.gz"
    )
    return dct


def getPairs_for_mapping(wildcards):
    dct = {}
    dct["Case"] = "mapStat/" + PAIRS[wildcards.sample][0] + ".flagstat"
    dct["Control"] = "mapStat/" + PAIRS[wildcards.sample][1] + ".flagstat"
    return dct


###############################################################
#####################主进程##################################
###############################################################


rule all:
    input:
        expand("extractReads/{sample}", sample=SAMPLES),
        expand("{sample}.html", sample=PAIRS.keys()),


rule ReadTrim:
    input:
        getFastqFiles,
    output:
        "cleanData/spCas9_{sample}_1_clean.fq.gz",
        "cleanData/spCas9_{sample}_2_clean.fq.gz",
    threads: config["threads"]
    singularity:
        "softwares/fastp-0.23.4.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "upstream_process"
    shell:
        "fastp -w {threads} -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]}"


rule CleanReadQC:
    input:
        "cleanData/spCas9_{sample}_1_clean.fq.gz",
        "cleanData/spCas9_{sample}_2_clean.fq.gz",
    output:
        directory("fastQC/{sample}"),
    threads: config["threads"]
    singularity:
        "softwares/fastqc-0.12.1.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "upstream_process"
    shell:
        "if [ ! -d fastQC/{wildcards.sample} ]; then mkdir -p fastQC/{wildcards.sample}; fi &&"
        "fastqc -t {threads} -o {output} {input}"


rule MergeByOverlap:
    input:
        "cleanData/spCas9_{sample}_1_clean.fq.gz",
        "cleanData/spCas9_{sample}_2_clean.fq.gz",
        "fastQC/{sample}",
    output:
        "flash/{sample}/{sample}.extendedFrags.fastq.gz",
    threads: config["threads"]
    singularity:
        "softwares/flash2-2.2.00.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "upstream_process"
    shell:
        "if [ ! -d flash/{wildcards.sample} ]; then mkdir -p flash/{wildcards.sample}; fi &&"
        "flash2 -t {threads} -m {config[MergeByOverlap][min_overlap]} -z -o {wildcards.sample} -d flash/{wildcards.sample} {input[0]} {input[1]}"


rule AlignToReference:
    input:
        "flash/{sample}/{sample}.extendedFrags.fastq.gz",
    output:
        temp("alignData/{sample}_raw.sam"),
    threads: config["threads"]
    singularity:
        "softwares/bwa-0.7.18.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "upstream_process"
    shell:
        "bwa mem -t {threads} {config[ref]} {input} > {output}"


rule AlignSortBam:
    input:
        "alignData/{sample}_raw.sam",
    output:
        "alignData/{sample}_sort.bam",
    threads: config["threads"]
    singularity:
        "softwares/samtools-1.19.2.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "upstream_process"
    shell:
        "samtools view -bF 4 {input} |"
        "samtools sort -@ {threads} -o {output} -"


rule mapStat:
    input:
        "alignData/{sample}_sort.bam",
    output:
        "mapStat/{sample}.flagstat",
    threads: config["threads"]
    singularity:
        "softwares/samtools-1.19.2.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "upstream_process"
    shell:
        "samtools flagstat -@ {threads} {input} > {output}"


rule IndexAlignmentFiles:
    input:
        "alignData/{sample}_sort.bam",
    output:
        "alignData/{sample}_sort.bam.bai",
    threads: config["threads"]
    singularity:
        "softwares/samtools-1.19.2.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "upstream_process"
    shell:
        "samtools index -@ {threads} {input}"


rule ExtractReads:
    input:
        bai="alignData/{sample}_sort.bam.bai",
        bam="alignData/{sample}_sort.bam",
        flagstat="mapStat/{sample}.flagstat",
    output:
        directory("extractReads/{sample}"),
    singularity:
        "softwares/python3-latest.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "upstream_process"
    script:
        "scripts/extractReads.py"


rule enrich_depeletion:
    input:
        unpack(getPairs),
    output:
        "basicStat/{sample}.enrich_depletion.tsv",
    singularity:
        "softwares/python3-latest.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "downstream_process"
    script:
        "scripts/enrich_depletion.py"


rule StatCleanQC:
    input:
        unpack(getPairs_for_qc),
    output:
        "basicStat/{sample}.cleanQC.tsv",
    singularity:
        "softwares/python3-latest.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "downstream_process"
    script:
        "scripts/StatCleanQC.py"


rule StatOverlap:
    input:
        unpack(getPairs_for_overlap),
    output:
        "basicStat/{sample}.StatOverlap.tsv",
    singularity:
        "softwares/python3-latest.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "downstream_process"
    script:
        "scripts/StatOverlap.py"


rule StatFlagstat:
    input:
        unpack(getPairs_for_mapping),
    output:
        "basicStat/{sample}.StatFlags.tsv",
    singularity:
        "softwares/python3-latest.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "downstream_process"
    script:
        "scripts/StatFlagstat.py"


rule RemoveSyntheticErrors:
    input:
        unpack(getPairs),
    output:
        directory("removeSyntheticErrors/{sample}"),
    singularity:
        "softwares/julia-latest.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "downstream_process"
    script:
        "scripts/removeSyntheticErrors.jl"


rule RealignToTrap:
    input:
        "removeSyntheticErrors/{sample}",
    output:
        directory("realignToTrap/{sample}"),
    singularity:
        "softwares/julia-latest.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "downstream_process"
    script:
        "scripts/realignToTrap.jl"


rule CreateIndelCountFile:
    input:
        "realignToTrap/{sample}",
    output:
        "basicStat/{sample}.indelCount.tsv",
    singularity:
        "softwares/julia-latest.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "downstream_process"
    script:
        "scripts/createIndelCountFile.jl"


rule visualization_whole:
    input:
        "basicStat/{sample}.indelCount.tsv",
    output:
        "visualizations/{sample}.indel_dist.png",
        "visualizations/{sample}.clean_reads_dist.png",
    singularity:
        "softwares/python3-latest.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "downstream_process"
    script:
        "scripts/Visualization2.py"


rule visualization_each:
    input:
        "realignToTrap/{sample}",
    output:
        directory("visualizations/{sample}"),
    singularity:
        "softwares/python3-latest.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "downstream_process"
    script:
        "scripts/Visualization.py"


rule createReport:
    input:
        indelCount="basicStat/{sample}.indelCount.tsv",
        qcStatfile="basicStat/{sample}.cleanQC.tsv",
        enrichfile="basicStat/{sample}.enrich_depletion.tsv",
        overlapStatfile="basicStat/{sample}.StatOverlap.tsv",
        mappingStatfile="basicStat/{sample}.StatFlags.tsv",
        indel_dist="visualizations/{sample}.indel_dist.png",
        cleanread_dist="visualizations/{sample}.clean_reads_dist.png",
        quiltdir="visualizations/{sample}",
        profilerdir="visualizations/{sample}",
    output:
        "{sample}.html",
    singularity:
        "softwares/python3-latest.sif"
    envmodules:
        "singularity/4.1.2",
    group:
        "downstream_process"
    script:
        "scripts/makeReport.py"
