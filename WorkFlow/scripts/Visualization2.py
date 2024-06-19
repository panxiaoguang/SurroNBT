import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def calculate_eff(file_path):
    df = pd.read_csv(file_path, sep="\t")
    df.columns = ["id", "clean", "indel"]
    df["eff"] = df["indel"] / df["clean"] * 100
    return df


def group_eff_and_plot(df, outfile):
    bins = [
        0,
        0.1,
        1,
        5,
        10,
        15,
        20,
        25,
        30,
        35,
        40,
        45,
        50,
        55,
        60,
        65,
        70,
        75,
        80,
        85,
        90,
        95,
        100,
    ]
    labels = [
        "0%-0.1%",
        "0.1%-1%",
        "1%-5%",
        "5%-10%",
        "10%-15%",
        "15%-20%",
        "20%-25%",
        "25%-30%",
        "30%-35%",
        "35%-40%",
        "40%-45%",
        "45%-50%",
        "50%-55%",
        "55%-60%",
        "60%-65%",
        "65%-70%",
        "70%-75%",
        "75%-80%",
        "80%-85%",
        "85%-90%",
        "90%-95%",
        "95%-100%",
    ]

    df["group"] = pd.cut(df["eff"], bins=bins, labels=labels, include_lowest=True)
    group_counts = df["group"].value_counts().sort_index()
    plt.figure(figsize=(6.6, 3.1))
    sns.barplot(
        x=group_counts.index, y=group_counts.values, color="#377eb8", edgecolor="black"
    )
    plt.xlabel("")
    plt.ylabel("Count")
    plt.title("")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)


def plot_log2_clean_distribution(file_path, output_file):
    df = pd.read_csv(file_path, sep="\t")
    df.columns = ["id", "clean", "indel"]
    plt.figure(figsize=(6.5, 2.8))
    sns.histplot(df, color="#377eb8", x="clean", log_scale=2)
    plt.xlabel("log2(clean)")
    plt.ylabel("Frequency")
    plt.title("Distribution of log2(clean)")
    plt.tight_layout()
    plt.savefig(output_file)


def main(input_file, indel_eff_out, clean_read_out):
    data = calculate_eff(input_file)
    group_eff_and_plot(data, indel_eff_out)
    plot_log2_clean_distribution(input_file, clean_read_out)


if __name__ == "__main__":
    main(snakemake.input[0], snakemake.output[0], snakemake.output[1])
