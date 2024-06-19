import gzip


def count_reads_in_fastq_gzip(file_path):
    line_count = 0
    with gzip.open(file_path, "rt") as f:
        for _ in f:
            line_count += 1
    reads_count = line_count // 4
    return reads_count


if __name__ == "__main__":
    reads_in_case = count_reads_in_fastq_gzip(snakemake.input["Case"])
    reads_in_control = count_reads_in_fastq_gzip(snakemake.input["Control"])
    with open(snakemake.output[0], "w") as f:
        f.write("reads_in_case\treads_in_control\n")
        f.write(f"{reads_in_case}\t{reads_in_control}\n")
