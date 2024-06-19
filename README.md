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
snakemake -s SurroNPROT.smk --cores 10 --use-singularity --singularity-args "--cleanenv --no-home"
```
You should prepare your input fasta files and put them in the directory: `WorkFlow/rawData/`
and modify the `WorkFlow/config/samples.tsv` according to the example file as well as the `WorkFlow/config/config.yaml` file.

We also provide an example input file in `rawData/` dir and the corresponding configuration file in `config/` dir. If your would like to try the example input, please delete all directories except `rawData/`, `config/ /references /scripts/`. Of course, the two files `SurroNPROT.smk` and `util.py` should be kept as well.

## Dependencies (Singularity container)
Please download all softwares for this pipeline in this link:
