from jinja2 import Environment, PackageLoader
import os
import pandas as pd
import glob


def number_of_passed_gRNAs(indelCount, libraryCount):
    if not os.path.exists(indelCount):
        raise FileNotFoundError(f"{indelCount} dosen't exists!")
    total_count = 0
    with open(indelCount, "r") as f:
        for line in f:
            if line.startswith("Label"):
                continue
            else:
                _, count_str, _ = line.strip().split("\t")
                count = int(count_str)
                if count >= 50:
                    total_count += 1
    return total_count, round(total_count / libraryCount * 100, 2)


def mean_indel_effs(indelCount):
    if not os.path.exists(indelCount):
        raise FileNotFoundError(f"{indelCount} dosen't exists!")
    df = pd.read_csv(indelCount, sep="\t")
    df.columns = ["gRNA", "totalCount", "indelCount"]
    ## calculate the mean indel efficiency: indelCount/totalCount
    df["indel_eff"] = df["indelCount"] / df["totalCount"]
    return round(df["indel_eff"].mean(), 2)


def get_quilts(inputdir):
    files = glob.glob(inputdir + "/quilt_plot_*.png")
    labels = [
        f.split("/")[-1].replace("quilt_plot_", "").replace(".png", "") for f in files
    ]
    return dict(zip(labels, files))


def get_profilers(inputdir):
    files = glob.glob(inputdir + "/alignment_plot_*.png")
    labels = [
        f.split("/")[-1].replace("alignment_plot_", "").replace(".png", "")
        for f in files
    ]
    return dict(zip(labels, files))


def main(
    indelCount,
    libraryCount,
    qcStatfile,
    enrichfile,
    overlapStatfile,
    mappingStatfile,
    quiltdir,
    profilerdir,
    indel_dist,
    cleanread_dist,
    outfile,
):
    env = Environment(loader=PackageLoader("report", "templates"))
    template = env.get_template("index.html")

    gRNA_passed, gRNA_passed_ratio = number_of_passed_gRNAs(indelCount, libraryCount)
    indeleff = mean_indel_effs(indelCount)
    df_enrich = pd.read_csv(enrichfile, sep="\t")
    df_qc = pd.read_csv(qcStatfile, sep="\t")
    df_flash = pd.read_csv(overlapStatfile, sep="\t")
    df_mapped = pd.read_csv(mappingStatfile, sep="\t")
    quilts = get_quilts(quiltdir)
    profilers = get_profilers(profilerdir)
    total_reads_case = max([int(x) for x in df_qc["totalSeq_case"][0].split("/")])
    total_reads_control = max([int(x) for x in df_qc["totalSeq_control"][0].split("/")])
    context = {
        "gRNA_passed": gRNA_passed,
        "gRNA_passed_ratio": gRNA_passed_ratio,
        "mean_indel_eff": indeleff,
        "enrich_depeltion": df_enrich["filter"].sum(),
        "Reads_Count_Case": df_qc["totalSeq_case"][0],
        "Reads_Count_Control": df_qc["totalSeq_control"][0],
        "Case_low": df_qc["failSeq_case"][0],
        "Control_low": df_qc["failSeq_control"][0],
        "mean_gc_Case": df_qc["gc_content_case"][0],
        "mean_gc_Control": df_qc["gc_content_control"][0],
        "Case_flash": df_flash["reads_in_case"][0] / total_reads_case * 100,
        "Control_flash": df_flash["reads_in_control"][0] / total_reads_control * 100,
        "Case_mapped": df_mapped["map_case"][0],
        "Control_mapped": df_mapped["map_control"][0],
        "indel_dist": indel_dist,
        "cleanread_dist": cleanread_dist,
        "quilts": quilts,
        "profiles": profilers,
    }
    with open(outfile, "w") as f:
        f.write(template.render(**context))


if __name__ == "__main__":
    indelCount = snakemake.input["indelCount"]
    libraryCount = snakemake.config["FinalReport"]["libraryCount"]
    qcStatfile = snakemake.input["qcStatfile"]
    enrichfile = snakemake.input["enrichfile"]
    overlapStatfile = snakemake.input["overlapStatfile"]
    mappingStatfile = snakemake.input["mappingStatfile"]
    quiltdir = snakemake.input["quiltdir"]
    profilerdir = snakemake.input["profilerdir"]
    indel_dist = snakemake.input["indel_dist"]
    cleanread_dist = snakemake.input["cleanread_dist"]
    outfile = snakemake.output[0]
    main(
        indelCount,
        libraryCount,
        qcStatfile,
        enrichfile,
        overlapStatfile,
        mappingStatfile,
        quiltdir,
        profilerdir,
        indel_dist,
        cleanread_dist,
        outfile,
    )
