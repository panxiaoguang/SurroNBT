configfile: "config/config.yaml"

import pandas as pd

def getSamples(iptFile):
    df = pd.read_table(iptFile)
    return df.SampleName.tolist()

def getPair(iptFile):
    df = pd.read_table(iptFile)
    dct = {}
    for row in df.to_dict(orient="records"):
        if isinstance(row["ControlName"],str):
            dct[row["SampleName"]] = [row["SampleName"], row["ControlName"]]
    return dct

SAMPLES = getSamples(config["samples"])
PAIRS = getPair(config["samples"])

def getFastqFiles(wildcards):
    df = pd.read_table(config["samples"])
    df = df[df.SampleName == wildcards.sample]
    if df.shape[0] != 1:
        raise Exception("Sample name not unique or not found in samples.tsv")
    fqs = [df.loc[1,"FastQ1"], df.loc[1,"FastQ2"]]
    return fqs

def getPairs(wildcards):
    dct = {}
    dct['Case'] = "extractReads/" + PAIRS[wildcards.sample][0]
    dct['Control'] = "extractReads/" + PAIRS[wildcards.sample][1]
    return dct

rule all:
    input:
        expand("extractReads/{sample}",sample=SAMPLES),
        expand("indelCounts/{sample}.tsv", sample = PAIRS.keys())

rule ReadTrim:
    input:
        getFastqFiles
    output:
        "cleanData/spCas9_{sample}_1_clean.fq.gz",
        "cleanData/spCas9_{sample}_2_clean.fq.gz"
    threads: config["threads"]
    shell:
        "fastp -w {threads} -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]}"

rule CleanReadQC:
    input:
        "cleanData/spCas9_{sample}_1_clean.fq.gz",
        "cleanData/spCas9_{sample}_2_clean.fq.gz"
    output:
        directory("fastQC/{sample}")
    threads: config["threads"]
    shell:
        "if [ ! -d fastQC/{wildcards.sample} ]; then mkdir -p fastQC/{wildcards.sample}; fi &&"
        "fastqc -t {threads} -o {output} {input}"

rule MergeByOverlap:
    input:
        "cleanData/spCas9_{sample}_1_clean.fq.gz",
        "cleanData/spCas9_{sample}_2_clean.fq.gz",
        "fastQC/{sample}"
    output:
        "flash/{sample}/{sample}.extendedFrags.fastq.gz"
    threads: config["threads"]
    shell:
        "if [ ! -d flash/{wildcards.sample} ]; then mkdir -p flash/{wildcards.sample}; fi &&"
        "flash -t {threads} -m {config[MergeByOverlap][min_overlap]} -z -o {wildcards.sample} -d flash/{wildcards.sample} {input[0]} {input[1]}"
rule AlignToReference:
    input:
        "flash/{sample}/{sample}.extendedFrags.fastq.gz"
    output:
        "alignData/{sample}.sort.bam"
    threads: config["threads"]
    shell:
        "bwa mem -t {threads} {config[ref]} {input} |"
        "samtools view -bF 4 - |"
        "samtools sort -@ {threads} -o {output}"

rule IndexAlignmentFiles:
    input:
        "alignData/{sample}.sort.bam"
    output:
        "alignData/{sample}.sort.bam.bai"
    threads: config["threads"]
    shell:
        "samtools index -@ {threads} {input}"

rule ExtractReads:
    input:
        bai = "alignData/{sample}.sort.bam.bai",
        bam = "alignData/{sample}.sort.bam"
    output:
        directory("extractReads/{sample}")
    script:
        "scripts/extractReads.py"

rule RemoveSyntheticErrors:
    input:
        unpack(getPairs)
    output:
        directory("removeSyntheticErrors/{sample}")
    script:
        "scripts/removeSyntheticErrors.jl"

rule RealignToTrap:
    input:
        "removeSyntheticErrors/{sample}"
    output:
        directory("realignToTrap/{sample}")
    script:
        "scripts/realignToTrap.jl"

rule CreateIndelCountFile:
    input:
        "realignToTrap/{sample}"
    output:
        "indelCounts/{sample}.tsv"
    script:
        "scripts/createIndelCountFile.jl"
