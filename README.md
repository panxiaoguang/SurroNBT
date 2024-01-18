## Welcome to SurroNPROT official repository!


## MainScripts

This is a directory that contains all scripts used in the SurroNPROT protocol.
please find them in the directory: `MainScripts/`   

### Scripts for library design

- `on-target_design.jl`: This script is used to design on-target SURRO_seq library for given regions of sgRNA.
- `utils/find_locations.jl`: This script is used to find the locations of given sgRNAs in the genome by alignment(bwa).
- `off-target_Artificial_design.jl`: This script is used to design off-target random mismatches library for given on-target sgRNA regions( get location by `find_locations.jl`).
- `off-target_genome_design.jl`: This script is used to design off-target SURRO_seq library for given on-target sgRNAs.
- `off-target_PAM_design.jl `: This script is used to design off-target random PAM only library for given on-target sgRNA regions( get location by `find_locations.jl`).
- `utils/remove_BsmBI.jl`: This script is used to remove BsmBI restriction enzyme sites in the designed library.

### Scripts for library analysis

- `get_seq.py `: This script is used to get the sequences(or we can call it library seperation) according to desgined library.
- `remove_synthesis_error.jl`: This script is used to remove synthesis error.
- `get_indel_profile.jl`: This script is used to get indel profile.
- `remove_sequencing_error.jl`: This script is used to remove sequencing error and get the final indel counts.


## One-key workflow

We also provide a one-key workflow for the whole SurroNPROT protocol. Please find it in the directory: `WorkFlow/`
It is a snakemake workflow, and you can run it by the following command:

```bash
cd WorkFlow/
snakemake -s SurroNPROT.smk --cores 10 
```

You should prepare your input fasta files and put them in the directory: `WorkFlow/rawData/`
and modify the `WorkFlow/config/samples.tsv` according to the example file as well as the `WorkFlow/config/config.yaml` file.

## Dependencies

- [Julia](https://julialang.org/) (>= 1.8.1)
- [Python](https://www.python.org/) (>= 3.6.0)
- [cas-offinder](https://github.com/snugel/cas-offinder/releases) (>= 2.4.1)
- [bwa](https://github.com/lh3/bwa) (>= 0.7.17)
- [samtools](https://github.com/samtools/samtools) (>= 1.13)
- [fastQC](https://github.com/s-andrews/FastQC)
- [fastp](https://github.com/OpenGene/fastp)
- [flash](https://ccb.jhu.edu/software/FLASH)
- [snakemake](https://snakemake.readthedocs.io/en/stable/)

### Julia packages
- BioSequences
- BioAlignments
- StatsBase
- FASTX
- Fire
- CSV
- DataFrames

### Python packages
- pysam
- numpy