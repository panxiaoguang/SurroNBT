import glob
import zipfile


def get_sequencing_stat(fastqc_path):
    total_reads = []
    reads_failed = []
    gc_content = []
    fastqc_zip_paths = glob.glob(fastqc_path + "/" + "*_fastqc.zip")
    for fastqc_zip_path in fastqc_zip_paths:
        with zipfile.ZipFile(fastqc_zip_path, "r") as zip_ref:
            # 找到 summary.txt 和 fastqc_data.txt 文件
            data_file = None
            for file_name in zip_ref.namelist():
                if file_name.endswith("fastqc_data.txt"):
                    data_file = file_name
            if data_file:
                with zip_ref.open(data_file) as f:
                    for line in f:
                        line = line.decode("utf-8")
                        if line.startswith("Total Sequences"):
                            total_reads.append(line.split("\t")[1].strip())
                        elif line.startswith("%GC"):
                            gc_content.append(line.split("\t")[1].strip())
                        elif line.startswith("Sequences flagged as poor quality"):
                            reads_failed.append(line.split("\t")[1].strip())
    return "/".join(total_reads), "/".join(reads_failed), "/".join(gc_content)


if __name__ == "__main__":
    fastqc_path_case = snakemake.input["Case"]
    fastqc_path_control = snakemake.input["Control"]
    totalSeq_case, failSeq_case, gc_content_case = get_sequencing_stat(fastqc_path_case)
    totalSeq_control, failSeq_control, gc_content_control = get_sequencing_stat(
        fastqc_path_control
    )
    with open(snakemake.output[0], "w") as f:
        f.write(
            "totalSeq_case\tfailSeq_case\tgc_content_case\ttotalSeq_control\tfailSeq_control\tgc_content_control\n"
        )
        f.write(
            f"{totalSeq_case}\t{failSeq_case}\t{gc_content_case}\t{totalSeq_control}\t{failSeq_control}\t{gc_content_control}\n"
        )
