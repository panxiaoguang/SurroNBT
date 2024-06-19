import plot_utils
import CRISPRessoPlot
import pandas as pd
import os
from pyfaidx import Fasta

customcolors = {
    "Substitution": "#0000FF",
    "Insertion": "#008000",
    "Deletion": "#FF0000",
    "A": "#7FC97F",
    "T": "#BEAED4",
    "C": "#FDC086",
    "G": "#FFFF99",
    "N": "#C8C8C8",
    "-": "#1E1E1E",
}


def create_quilt_plot(
    refseq, profiler, customcolors, output_prefiex, sgS, sgE, lftS, lftE, rgtS, rgtE
):
    editing_dict, indel_dict, total_reads = plot_utils.substitute_position(
        profiler, refseq, lftS, lftE, rgtS, rgtE
    )
    nuc_df = pd.DataFrame(editing_dict)
    nuc_df.columns = list(refseq)
    nuc_df = (
        nuc_df.divide(total_reads).reset_index().rename(columns={"index": "Nucleotide"})
    )
    nuc_df.insert(0, "Batch", "sgRNA")

    modiy_df = pd.DataFrame(indel_dict)
    modiy_df.columns = list(refseq)
    modiy_df = (
        modiy_df.divide(total_reads)
        .reset_index()
        .rename(columns={"index": "Modification"})
    )
    modiy_df.insert(0, "Batch", "sgRNA")
    CRISPRessoPlot.plot_nucleotide_quilt(
        nuc_df,
        modiy_df,
        output_prefiex,
        customcolors,
        sgRNA_intervals=[[sgS - 1, sgE - 1]],
        save_also_png=True,
    )


def create_alignment_heatmap(
    refseq,
    profiler,
    customcolors,
    output_prefiex,
    sgS,
    sgE,
    lftS,
    lftE,
    rgtS,
    rgtE,
    trapL,
):
    freq_table = plot_utils.make_frequency_table(
        profiler, lftS, lftE, rgtS, rgtE, trapL
    )
    freq_table = freq_table.reset_index().set_index("Aligned_Sequence")
    CRISPRessoPlot.plot_alleles_table(
        refseq,
        freq_table,
        output_prefiex,
        customcolors,
        sgRNA_intervals=[[sgS - 1, sgE - 1]],
        cutting_point=sgE - 3,
        SAVE_ALSO_PNG=True,
    )


def main(
    reference,
    inpath,
    outpath,
    leftE,
    sgS,
    sgE,
    lftS,
    lftE,
    rgtS,
    rgtE,
    customcolors,
    trapL,
):
    ref_reader = Fasta(reference)
    for record in ref_reader:
        ref_sequence = str(record)
        refseq = ref_sequence[leftE : (leftE + trapL)]
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        outprefix_quilt = os.path.join(outpath, f"quilt_plot_sgRNA_{record.name}")
        outprefix_alignment = os.path.join(
            outpath, f"alignment_plot_sgRNA_{record.name}"
        )
        in_file = os.path.join(inpath, f"{record.name}.indelProfile.tsv")
        create_quilt_plot(
            refseq,
            in_file,
            customcolors,
            outprefix_quilt,
            sgS,
            sgE,
            lftS,
            lftE,
            rgtS,
            rgtE,
        )
        create_alignment_heatmap(
            refseq,
            in_file,
            customcolors,
            outprefix_alignment,
            sgS,
            sgE,
            lftS,
            lftE,
            rgtS,
            rgtE,
            trapL,
        )


if __name__ == "__main__":
    main(
        snakemake.config["ref"],
        snakemake.input[0],
        snakemake.output[0],
        snakemake.config["ExtractReads"]["leftEnd"],
        snakemake.config["Visualization"]["sgRNA_start"],
        snakemake.config["Visualization"]["sgRNA_end"],
        snakemake.config["CreateIndelCountFile"]["leftStart"],
        snakemake.config["CreateIndelCountFile"]["leftEnd"],
        snakemake.config["CreateIndelCountFile"]["rightStart"],
        snakemake.config["CreateIndelCountFile"]["rightEnd"],
        customcolors,
        snakemake.config["ExtractReads"]["trapLength"],
    )
