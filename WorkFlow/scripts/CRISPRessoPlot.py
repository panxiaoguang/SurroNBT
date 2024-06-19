import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
import numpy as np
from collections import defaultdict
import re
from matplotlib import colors as colors_mpl
import matplotlib.gridspec as gridspec
import seaborn as sns


def setMatplotlibDefaults():
    font = {'size': 22}
    matplotlib.rc('font', **font)
    matplotlib.use('AGG')
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    matplotlib.rcParams["font.sans-serif"] = ["Arial",
                                              "Liberation Sans", "Bitstream Vera Sans"]
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams['axes.facecolor'] = 'white'
    sns.set(style='white', font_scale=2.2)
    plt.ioff()


setMatplotlibDefaults()


def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def get_nuc_color(nuc, alpha):
    def get_color(x, y, z): return (x/255.0, y/255.0, z/255.0, alpha)
    if nuc == "A":
        return get_color(127, 201, 127)
    elif nuc == "T":
        return get_color(190, 174, 212)
    elif nuc == "C":
        return get_color(253, 192, 134)
    elif nuc == "G":
        return get_color(255, 255, 153)
    elif nuc == "N":
        return get_color(200, 200, 200)
    elif nuc == "INS":
        #        return get_color(185,219,253)
        #        return get_color(177,125,76)
        return get_color(193, 129, 114)
    elif nuc == "DEL":
        # return get_color(177,125,76)
        #        return get_color(202,109,87)
        return get_color(193, 129, 114)
    elif nuc == "-":
        # return get_color(177,125,76)
        #        return get_color(202,109,87)
        return get_color(30, 30, 30)
    else:  # return a random color (that is based on the nucleotide given)
        charSum = 0
        for char in nuc.upper():
            thisval = ord(char) - 65  # 'A' is 65
            if thisval < 0 or thisval > 90:
                thisval = 0
            charSum += thisval
        charSum = (charSum/len(nuc))/90.0

        return (charSum, (1-charSum), (2*charSum*(1-charSum)))


def get_color_lookup(nucs, alpha, custom_colors=None):
    if custom_colors is None:
        colorLookup = {}
        for nuc in nucs:
            colorLookup[nuc] = get_nuc_color(nuc, alpha)
        return colorLookup
    else:
        def get_color(x, y, z): return (x / 255.0, y / 255.0, z / 255.0, alpha)
        colors = {}
        for nuc in nucs:
            if nuc == 'INS':
                rgb = (193, 129, 114)
            else:
                rgb = hex_to_rgb(custom_colors[nuc])
            colors[nuc] = get_color(rgb[0], rgb[1], rgb[2])
        return colors


def get_rows_for_sgRNA_annotation(sgRNA_intervals, amp_len):
    """
    Returns an array specifying the row number that an sgRNA should be plotted on in order to avoid overlap
    params:
    sgRNA_intervals: array of x coordinate tuples of start and stop
    amp_len: length of amplicon
    returns: sgRNA_plot_rows: list of index on which row to plot

    """
    # figure out how many rows are needed to show all sgRNAs
    # which row each sgRNA should be plotted on
    sgRNA_plot_rows = [0]*len(sgRNA_intervals)
    sgRNA_plot_occupancy = []  # which idxs are already filled on each row
    sgRNA_plot_occupancy.append([])
    for idx, sgRNA_int in enumerate(sgRNA_intervals):
        this_sgRNA_start = max(0, sgRNA_int[0])
        this_sgRNA_end = min(sgRNA_int[1], amp_len - 1)
        curr_row = 0
        if this_sgRNA_start > amp_len or this_sgRNA_end < 0:
            sgRNA_plot_rows[idx] = curr_row
            continue
        while len(np.intersect1d(sgRNA_plot_occupancy[curr_row], range(this_sgRNA_start, this_sgRNA_end))) > 0:
            next_row = curr_row + 1
            if not next_row in sgRNA_plot_occupancy:
                sgRNA_plot_occupancy.append([])
            curr_row = next_row
        sgRNA_plot_rows[idx] = curr_row
        sgRNA_plot_occupancy[curr_row].extend(
            range(this_sgRNA_start, this_sgRNA_end))
    return (np.subtract(max(sgRNA_plot_rows), sgRNA_plot_rows))


def add_sgRNA_to_ax(ax, sgRNA_intervals, sgRNA_y_start, sgRNA_y_height, amp_len, x_offset=0, sgRNA_mismatches=None, sgRNA_names=None, sgRNA_rows=None, font_size=None, clip_on=True, label_at_zero=False):
    """
    Adds sgRNA to plot ax
    params:
    ax: ax to add sgRNA to
    sgRNA_intervals: array of x coordinate tuples of start and stop
    sgRNA_y_start: y coordinate of sgRNA(s)
    sgRNA_y_height: y height of sgRNA(s)
    amp_len: length of amplicon
    x_offset: amount to move sgRNAs in x direction -- if labels or annotations are to the left of the plot, this may have to be non-zero. e.g. if the reference starts at x=2, set this to 2
    sgRNA_mismatches: array (for each sgRNA_interval) of locations in sgRNA where there are mismatches
    sgRNA_names: array (for each sgRNA_interval) of names of sgRNAs (otherwise empty)
    sgRNA_rows: which row to plot the sgRNA on so they don't overlap (y-axis-wise)
    clip_on: matplotlib parameter for whether sgRNAs should be drawn outside of clipping bounds (if sgRNAs aren't showing up, try setting this to False)
    label_at_zero: whether first sgRNA should be forced to be at 0 instead of off the plot to the left beyond 0 (some plots are ok with this)
    """

    if font_size is None:
        font_size = matplotlib.rcParams['font.size']

    # figure out how many rows are needed to show all sgRNAs
    if sgRNA_rows is None:
        sgRNA_rows = [0]*len(sgRNA_intervals)
    max_sgRNA_row = max(sgRNA_rows)+1
    this_sgRNA_y_height = sgRNA_y_height/float(max_sgRNA_row)

    min_sgRNA_x = None  # keep track of left-most sgRNA
    # whether to label left-most sgRNA (set to false if label another sgRNA (e.g. with sgRNA_name))
    label_left_sgRNA = True
    for idx, sgRNA_int in enumerate(sgRNA_intervals):
        this_sgRNA_start = max(0, sgRNA_int[0])
        this_sgRNA_end = min(sgRNA_int[1], amp_len - 1)
        if this_sgRNA_start > amp_len or this_sgRNA_end < 0:
            continue
        this_sgRNA_y_start = sgRNA_y_start + \
            this_sgRNA_y_height*sgRNA_rows[idx]
        ax.add_patch(
            patches.Rectangle((x_offset+this_sgRNA_start, this_sgRNA_y_start), 1+this_sgRNA_end -
                              this_sgRNA_start, this_sgRNA_y_height, facecolor=(0, 0, 0, 0.15), clip_on=clip_on)
        )

        # if plot has trimmed the sgRNA, add a mark
        if (this_sgRNA_start) != sgRNA_int[0]:
            ax.add_patch(
                # patches.Rectangle((x_offset + 0.1+this_sgRNA_start, sgRNA_y_start), 0.1, sgRNA_y_height,facecolor='w',clip_on=clip_on)
                patches.Rectangle((x_offset + 0.1+this_sgRNA_start, this_sgRNA_y_start),
                                  0.1, this_sgRNA_y_height, facecolor='w', clip_on=clip_on)
            )
        if this_sgRNA_end != sgRNA_int[1]:
            ax.add_patch(
                patches.Rectangle((x_offset + 0.8+this_sgRNA_end, this_sgRNA_y_start),
                                  0.1, this_sgRNA_y_height, facecolor='w', clip_on=clip_on)
            )

        if sgRNA_mismatches is not None:
            this_sgRNA_mismatches = sgRNA_mismatches[idx]
            for mismatch in this_sgRNA_mismatches:
                mismatch_plot_pos = sgRNA_int[0] + mismatch
                if mismatch_plot_pos > 0 and mismatch_plot_pos < amp_len - 1:
                    ax.add_patch(
                        patches.Rectangle((x_offset + mismatch_plot_pos, this_sgRNA_y_start),
                                          1, this_sgRNA_y_height, facecolor='r', clip_on=clip_on)
                    )

        # set left-most sgrna start
        if not min_sgRNA_x:
            min_sgRNA_x = this_sgRNA_start
        if this_sgRNA_start < min_sgRNA_x:
            min_sgRNA_x = this_sgRNA_start
        if sgRNA_names is not None and idx < len(sgRNA_names) and sgRNA_names[idx] != "":
            if (label_at_zero and x_offset + this_sgRNA_start < len(sgRNA_names[idx])*0.66):
                ax.text(0, this_sgRNA_y_start + this_sgRNA_y_height/2,
                        sgRNA_names[idx] + " ", horizontalalignment='left', verticalalignment='center', fontsize=font_size)
            else:
                ax.text(x_offset+this_sgRNA_start, this_sgRNA_y_start + this_sgRNA_y_height/2,
                        sgRNA_names[idx] + " ", horizontalalignment='right', verticalalignment='center', fontsize=font_size)
            label_left_sgRNA = False  # already labeled at least one sgRNA

    if min_sgRNA_x is not None and label_left_sgRNA:
        if (label_at_zero and x_offset + min_sgRNA_x < 5):
            ax.text(0, this_sgRNA_y_start + this_sgRNA_y_height/2, 'sgRNA ',
                    horizontalalignment='left', verticalalignment='center', fontsize=font_size)
        else:
            ax.text(x_offset+min_sgRNA_x, this_sgRNA_y_start + this_sgRNA_y_height/2, 'sgRNA ',
                    horizontalalignment='right', verticalalignment='center', fontsize=font_size)


def plot_nucleotide_quilt(nuc_pct_df, mod_pct_df, fig_filename_root, custom_colors, save_also_png=False, sgRNA_intervals=None, min_text_pct=0.5, max_text_pct=0.95, quantification_window_idxs=None, sgRNA_names=None, sgRNA_mismatches=None, shade_unchanged=True, group_column='Batch', **kwargs):
    """
    Plots a nucleotide quilt with each square showing the percentage of each base at that position in the reference
    nuc_pct_df: dataframe with percents of each base (ACTGN-) at each position
    mod_pct_df: dataframe with percents of modifications at each position (this function uses 'Insertions_Left' to plot insertions)
    fig_filename_root: filename root (will add .pdf or .png)
    save_also_png: whether png should also be saved
    sgRNA_intervals: ranges for sgRNA annotation on plot
    sgRNA_names: names to annotate sgRNAs with (if None, will just label left sgRNA with 'sgRNA')
    sgRNA_mismatches: locations in the sgRNA where there are mismatches from an original guide (flexiguides)
    quantification_window_idxs: indices for quantification window annotation on plot
    min_text_pct: add text annotation if the percent is greater than this number
    max_text_pct: add text annotation if the percent is less than this number
    shade_unchanged: if true, unchanged/reference nucleotides will be shaded (only changes with regard to reference will be dark)
    group_column: If multiple samples are given, they are grouped by this column
    """
    plotPct = 0.9  # percent of vertical space to plot in (the rest will be white)
    # if value is less than this, it won't plot the rectangle (with white boundary)
    min_plot_pct = 0.01

    if float(nuc_pct_df.iloc[1, 2]) > 1:
        raise Exception(
            'Expecting nucleotide percentage. Instead, got numbers in nuc_pct_df: ' + str(nuc_pct_df.iloc[1, 2]))

    nrows = nuc_pct_df.shape[0]
    # Batch, Nucleotide, nuc1, nuc2, nuc3 ...
    amp_len = nuc_pct_df.shape[1] - 2
    nucs = nuc_pct_df.Nucleotide.unique()
    nNucs = len(nucs)
    nSamples = int(nrows / nNucs)
    samplesList = []
    for i in range(nSamples):  # iterate over all samples
        sample_row_start = nNucs * i
        samplesList.append(nuc_pct_df.iloc[sample_row_start, 0])

    # make a color map of fixed colors
    color_lookup = get_color_lookup(
        ['A', 'T', 'C', 'G', 'N', 'INS', '-'], alpha=1, custom_colors=custom_colors)
    unchanged_color_lookup = get_color_lookup(['A', 'T', 'C', 'G', 'N', 'INS', '-'], alpha=0.3,
                                              custom_colors=custom_colors)

    # fig = plt.figure(figsize=(amp_len/2.0,nSamples*2))
    # fig = plt.figure(figsize=(amp_len,nSamples))
    # fig = plt.figure(figsize=(amp_len,nSamples*2))
    # fig = plt.figure(figsize=(amp_len,(nSamples+1)*2))
    fig, ax = plt.subplots(figsize=((amp_len+10)/2.0, (nSamples+1)*2))

    # remove box around plot
    for spine in ax.spines.values():
        spine.set_visible(False)

    if not shade_unchanged:  # shade all nucs equally
        # iterate over all nucleotide positions in the sequence (0=Batch, 1=Nucleotide, so start at 2)
        for pos_ind in range(2, amp_len+2):
            x_start = pos_ind
            x_end = pos_ind + 1
            for i in range(nSamples):  # iterate over all samples
                sample_row_start = nNucs * i
                y_start = nSamples - i
                sumPct = 0
                # iterate over each nucleotide at this position in this sample
                for nuc_ind in range(nNucs):
                    pct = float(
                        nuc_pct_df.iloc[sample_row_start + nuc_ind, pos_ind])
                    sumPct += pct
                    if pct > min_plot_pct:
                        obs_pct = pct * plotPct
                        curr_nuc = nuc_pct_df.iloc[sample_row_start + nuc_ind, 1]
                        ax.add_patch(
                            patches.Rectangle(
                                (x_start, y_start), x_end-x_start, obs_pct, facecolor=color_lookup[curr_nuc], edgecolor='w')
                        )
                        if pct > min_text_pct and pct < max_text_pct:
                            ax.text(x_start+0.55, y_start + obs_pct/2.0, format(pct*100, '.1f'),
                                    horizontalalignment='center', verticalalignment='center', rotation=90)
                        y_start += obs_pct

    else:  # shade unchanged bases
        ref_seq = nuc_pct_df.columns.values
        # iterate over all nucleotide positions in the sequence (0=Batch, 1=Nucleotide, so start at 2)
        for pos_ind in range(2, amp_len+2):
            x_start = pos_ind
            x_end = pos_ind + 1
            for i in range(nSamples):  # iterate over all samples
                sample_row_start = nNucs * i
                y_start = nSamples - i
                sumPct = 0
                # iterate over each nucleotide at this position in this sample
                for nuc_ind in range(nNucs):
                    pct = float(
                        nuc_pct_df.iloc[sample_row_start + nuc_ind, pos_ind])
                    sumPct += pct
                    if pct > min_plot_pct:
                        obs_pct = pct * plotPct
                        curr_nuc = nuc_pct_df.iloc[sample_row_start + nuc_ind, 1]
                        if curr_nuc == ref_seq[pos_ind]:  # if is reference
                            ax.add_patch(
                                patches.Rectangle((x_start, y_start), x_end-x_start, obs_pct,
                                                  facecolor=unchanged_color_lookup[curr_nuc], edgecolor='w')
                            )
                        else:
                            ax.add_patch(
                                patches.Rectangle(
                                    (x_start, y_start), x_end-x_start, obs_pct, facecolor=color_lookup[curr_nuc], edgecolor='w')
                            )

                        if pct > min_text_pct and pct < max_text_pct:
                            ax.text(x_start+0.55, y_start + obs_pct/2.0, format(pct*100, '.1f'),
                                    horizontalalignment='center', verticalalignment='center', rotation=90)
                        y_start += obs_pct

    mod_pct_df_indexed = mod_pct_df.set_index([group_column, 'Modification'])
    # add insertions
    # iterate over all nucleotide positions in the sequence (0=Batch, 1=Modification, so start at 2)
    for pos_ind in range(2, amp_len+1):
        x_start = pos_ind + 0.7
        x_end = pos_ind + 1.3
        for i in range(nSamples):  # iterate over all samples
            sampleName = samplesList[i]

            sample_row_start = nNucs * i
            y_start = nSamples - i

            ins_pct = float(
                mod_pct_df_indexed.loc[sampleName, 'Insertions_Left'].iloc[pos_ind-2])

            if ins_pct > min_plot_pct:
                obs_pct = ins_pct * plotPct
                ax.add_patch(
                    patches.Rectangle((x_start, y_start), x_end-x_start,
                                      obs_pct, facecolor=color_lookup['INS'], edgecolor='w')
                )
                if ins_pct > min_text_pct and ins_pct < max_text_pct:
                    ax.text(x_start+0.15, y_start + obs_pct/2.0, format(ins_pct*100, '.1f'),
                            horizontalalignment='center', verticalalignment='center', rotation=90)

    # draw black box around each sample
    for i in range(nSamples):
        y_start = nSamples - i
        ax.add_patch(
            patches.Rectangle((2, y_start), amp_len, plotPct,
                              facecolor='none', edgecolor='black')
        )

    # draw reference sequence
    ref_y_start = 0.5
    ref_y_height = 0.4
    ref_seq = nuc_pct_df.columns.values
    # iterate over all nucleotide positions in the sequence (0=Batch, 1=Nucleotide, so start at 2)
    for pos_ind in range(2, amp_len+2):
        ax.add_patch(
            patches.Rectangle((pos_ind, ref_y_start), 1, ref_y_height,
                              facecolor=color_lookup[ref_seq[pos_ind]], edgecolor='w')
        )
        ax.text(pos_ind+0.5, ref_y_start + ref_y_height/2.3,
                ref_seq[pos_ind], horizontalalignment='center', verticalalignment='center')

    ax.tick_params(top=False, bottom=False, left=False,
                   right=False, labelleft=True, labelbottom=False)

    ax.set_yticks([ref_y_start + ref_y_height/2.0] +
                  [x+0.5 for x in range(1, nSamples+1)])
#    sampleLabs = list(nuc_pct_df.iloc[[((nSamples-1)-x)*nNucs for x in range(0,nSamples)],0]))
#    print(mod_pct_df)
#    sampleReadCounts = list(nuc_pct_df.iloc[[((nSamples-1)-x)*nNucs for x in range(0,nSamples)],0]))
    ax.set_yticklabels(['Reference'] + list(nuc_pct_df.iloc[[((nSamples-1)-x)
                       * nNucs for x in range(0, nSamples)], 0]), va='center')

    plot_y_start = ref_y_start - 0.1

    if sgRNA_intervals:
        sgRNA_rows = get_rows_for_sgRNA_annotation(sgRNA_intervals, amp_len)
        num_sgRNA_rows = max(sgRNA_rows) + 1
        sgRNA_y_height = num_sgRNA_rows * 0.3
        plot_y_start = ref_y_start - (sgRNA_y_height + 0.1)
        add_sgRNA_to_ax(ax, sgRNA_intervals, sgRNA_y_start=plot_y_start + 0.1, sgRNA_y_height=sgRNA_y_height-0.1,
                        amp_len=amp_len, x_offset=2, sgRNA_mismatches=sgRNA_mismatches, sgRNA_names=sgRNA_names, sgRNA_rows=sgRNA_rows)

    if quantification_window_idxs is not None and len(quantification_window_idxs) > 0:
        q_win_y_start = plot_y_start
        q_win_y_height = nSamples+1 - q_win_y_start

        q_list = sorted(list(quantification_window_idxs))

        lastStart = q_list[0]
        lastIdx = q_list[0]
        for idx in range(1, len(q_list)):
            if q_list[idx] == lastIdx + 1:
                lastIdx = q_list[idx]
            else:
                ax.add_patch(
                    patches.Rectangle((2+lastStart, q_win_y_start), 1+(lastIdx-lastStart), q_win_y_height,
                                      fill=None, edgecolor=(0, 0, 0, 0.25), linestyle=(0, (5, 2)), linewidth=2)
                )
                lastStart = q_list[idx]
                lastIdx = q_list[idx]
        ax.add_patch(
            patches.Rectangle((2+lastStart, q_win_y_start), 1+(lastIdx-lastStart), q_win_y_height,
                              fill=None, edgecolor=(0, 0, 0, 0.25), linestyle=(0, (5, 2)), linewidth=2)
        )

    ax.set_xlim([2, amp_len+3])
    ax.set_ylim([plot_y_start, nSamples+1.2])

    legend_patches = []
    for nuc in nucs:
        if nuc == "-":
            continue
        if nuc == "N":
            continue
        patch = patches.Patch(color=color_lookup[nuc], label=nuc)
        legend_patches.append(patch)

    n_tab = nuc_pct_df.loc[nuc_pct_df.Nucleotide == 'N']
    n_sum = n_tab.iloc[:, 2:len(n_tab.columns)].sum().sum()
    ins_tab = mod_pct_df.loc[mod_pct_df.Modification == 'Insertions_Left']
    ins_sum = ins_tab.iloc[:, 2:len(ins_tab.columns)].sum().sum()
    del_tab = mod_pct_df.loc[mod_pct_df.Modification == 'Deletions']
    del_sum = del_tab.iloc[:, 2:len(del_tab.columns)].sum().sum()
    if n_sum > 0:
        n_patch = patches.Patch(color=color_lookup['N'], label='N')
        legend_patches.append(n_patch)
    if ins_sum > 0:
        ins_patch = patches.Patch(color=color_lookup['INS'], label='Insertion')
        legend_patches.append(ins_patch)
    if del_sum > 0:
        del_patch = patches.Patch(color=color_lookup['-'], label='Deletion')
        legend_patches.append(del_patch)
    if quantification_window_idxs is not None and len(quantification_window_idxs) > 0:
        q_win_patch = patches.Patch(fill=None, edgecolor=(0, 0, 0, 0.25), linestyle=(
            0, (5, 2)), linewidth=2, label='Quantification window')
        legend_patches.append(q_win_patch)

    ax.legend(handles=legend_patches, loc='center left',
              ncol=1, bbox_to_anchor=(1, 0.5))

    # todo -- if the plot_around_cut is really small (e.g. 2) the plots are blown out of proportion.. this could be fixed here, but not easily
#    bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#    width = bbox.width
#    height = bbox.height
#    print('width is ' + str(width) + ' and height is ' + str(height))
#    if (width < 50):
#        print('setting here!!')
# fig.set_figwidth(50)
# fig.tight_layout(w_pad=0.5)
#        fig.tight_layout(h_pad=0.5)
#    else:
#        fig.tight_layout()

    fig.savefig(fig_filename_root+'.pdf', bbox_inches='tight')
    if save_also_png:
        fig.savefig(fig_filename_root+'.png',
                    bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)


def prep_alleles_table(df_alleles, reference_seq, MAX_N_ROWS, MIN_FREQUENCY):
    """
    Prepares a df of alleles for Plotting
    input:
    -df_alleles: pandas dataframe of alleles to plot
    -reference_seq: sequence of unmodified reference
    -MAX_N_ROWS: max number of rows to plot
    -MIN_FREQUENCY: min frequency for a row to be plotted
    returns:
    -X: list of numbers representing nucleotides of the allele
    -annot: list of nucleotides (letters) of the allele
    -y_labels: list of labels for each row/allele
    -insertion_dict: locations of insertions -- red squares will be drawn around these
    -per_element_annot_kws: annotations for each cell (e.g. bold for substitutions, etc.)
    -is_reference: list of booleans for whether the read is equal to the reference
    """
    dna_to_numbers = {'-': 0, 'A': 1, 'T': 2, 'C': 3, 'G': 4, 'N': 5}
    def seq_to_numbers(seq): return [dna_to_numbers[x] for x in seq]

    X = []
    annot = []
    y_labels = []
    insertion_dict = defaultdict(list)
    per_element_annot_kws = []
    is_reference = []

    re_find_indels = re.compile("(-*-)")
    idx_row = 0
    for idx, row in df_alleles[df_alleles['%Reads'] >= MIN_FREQUENCY][:MAX_N_ROWS].iterrows():
        X.append(seq_to_numbers(idx.upper()))
        annot.append(list(idx))

        has_indels = False
        for p in re_find_indels.finditer(row['Reference_Sequence']):
            has_indels = True
            insertion_dict[idx_row].append((p.start(), p.end()))

        y_labels.append('%.2f%% (%d reads)' % (row['%Reads'], row['#Reads']))
        if idx == reference_seq and not has_indels:
            is_reference.append(True)
        else:
            is_reference.append(False)

        idx_row += 1

        idxs_sub = [i_sub for i_sub in range(len(idx)) if
                    (row['Reference_Sequence'][i_sub] != idx[i_sub]) and
                    (row['Reference_Sequence'][i_sub] != '-') and
                    (idx[i_sub] != '-')]
        to_append = np.array([{}]*len(idx), dtype=object)
        to_append[idxs_sub] = {'weight': 'bold', 'color': 'black', 'size': 16}
        per_element_annot_kws.append(to_append)

    return X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference


class Custom_HeatMapper(sns.matrix._HeatMapper):

    def __init__(self, data, vmin, vmax, cmap, center, robust, annot, fmt,
                 annot_kws, per_element_annot_kws, cbar, cbar_kws,
                 xticklabels=True, yticklabels=True, mask=None):

        super(Custom_HeatMapper, self).__init__(data, vmin, vmax, cmap, center, robust, annot, fmt,
                                                annot_kws, cbar, cbar_kws,
                                                xticklabels, yticklabels, mask)

        if annot is not None:
            if per_element_annot_kws is None:
                self.per_element_annot_kws = np.empty_like(annot, dtype=object)
                self.per_element_annot_kws[:] = dict()
            else:
                self.per_element_annot_kws = per_element_annot_kws

    # add per element dict to syle the annotatiin
    def _annotate_heatmap(self, ax, mesh):
        """Add textual labels with the value in each cell."""
        mesh.update_scalarmappable()
        xpos, ypos = np.meshgrid(ax.get_xticks(), ax.get_yticks())

        for x, y, m, color, val, per_element_dict in zip(xpos.flat, ypos.flat,
                                                         mesh.get_array().flat, mesh.get_facecolors(),
                                                         self.annot_data.flat, self.per_element_annot_kws.flat):
            # print per_element_dict
            if m is not np.ma.masked:
                l = sns.utils.relative_luminance(color)
                text_color = ".15" if l > .408 else "w"
                annotation = ("{:" + self.fmt + "}").format(str(val))
                text_kwargs = dict(color=text_color, ha="center", va="center")
                text_kwargs.update(self.annot_kws)
                text_kwargs.update(per_element_dict)

                ax.text(x, y, annotation, **text_kwargs)

    # removed the colobar

    def plot(self, ax, cax, kws):
        """Draw the heatmap on the provided Axes."""
        # Remove all the Axes spines
        sns.utils.despine(ax=ax, left=True, bottom=True)

        # Draw the heatmap
        mesh = ax.pcolormesh(self.plot_data, vmin=self.vmin, vmax=self.vmax,
                             cmap=self.cmap, **kws)

        # Set the axis limits
        ax.set(xlim=(0, self.data.shape[1]), ylim=(0, self.data.shape[0]))

        # Add row and column labels
        ax.set(xticks=self.xticks, yticks=self.yticks)
        xtl = ax.set_xticklabels(self.xticklabels)
        ytl = ax.set_yticklabels(
            self.yticklabels, rotation="vertical", va='center')

        # Possibly rotate them if they overlap
        plt.draw()
        if sns.utils.axis_ticklabels_overlap(xtl):
            plt.setp(xtl, rotation="vertical")
        if sns.utils.axis_ticklabels_overlap(ytl):
            plt.setp(ytl, rotation="horizontal")

        # Add the axis labels
        ax.set(xlabel=self.xlabel, ylabel=self.ylabel)

        # Annotate the cells with the formatted values
        if self.annot:
            self._annotate_heatmap(ax, mesh)


def custom_heatmap(data, vmin=None, vmax=None, cmap=None, center=None, robust=False,
                   annot=None, fmt=".2g", annot_kws=None, per_element_annot_kws=None,
                   linewidths=0, linecolor="white",
                   cbar=True, cbar_kws=None, cbar_ax=None,
                   square=False, ax=None, xticklabels=True, yticklabels=True,
                   mask=None,
                   **kwargs):

    # Initialize the plotter object
    plotter = Custom_HeatMapper(data, vmin, vmax, cmap, center, robust, annot, fmt,
                                annot_kws, per_element_annot_kws, cbar, cbar_kws, xticklabels,
                                yticklabels, mask)

    # Add the pcolormesh kwargs here
    kwargs["linewidths"] = linewidths
    kwargs["edgecolor"] = linecolor

    # Draw the plot and return the Axes
    if ax is None:
        ax = plt.gca()
    if square:
        ax.set_aspect("equal")
    plotter.plot(ax, cbar_ax, kwargs)
    return ax


def plot_alleles_heatmap(
        reference_seq,
        fig_filename_root,
        X,
        annot,
        y_labels,
        insertion_dict,
        per_element_annot_kws,
        custom_colors,
        SAVE_ALSO_PNG=False,
        plot_cut_point=True,
        sgRNA_intervals=None,
        sgRNA_names=None,
        sgRNA_mismatches=None,
        cutting_point = 26.5,
        **kwargs):
    """
    Plots alleles in a heatmap (nucleotides color-coded for easy visualization)
    input:
    -reference_seq: sequence of reference allele to plot
    -fig_filename: figure filename to plot (not including '.pdf' or '.png')
    -X: list of numbers representing nucleotides of the allele
    -annot: list of nucleotides (letters) of the allele
    -y_labels: list of labels for each row/allele
    -insertion_dict: locations of insertions -- red squares will be drawn around these
    -per_element_annot_kws: annotations for each cell (e.g. bold for substitutions, etc.)
    -SAVE_ALSO_PNG: whether to write png file as well
    -plot_cut_point: if false, won't draw 'predicted cleavage' line
    -sgRNA_intervals: locations where sgRNA is located
    -sgRNA_mismatches: array (for each sgRNA_interval) of locations in sgRNA where there are mismatches
    -sgRNA_names: array (for each sgRNA_interval) of names of sgRNAs (otherwise empty)
    -custom_colors: dict of colors to plot (e.g. colors['A'] = (1,0,0,0.4) # red,blue,green,alpha )
    """
    plot_nuc_len = len(reference_seq)

    # make a color map of fixed colors
    alpha = 0.4
    A_color = get_nuc_color('A', alpha)
    T_color = get_nuc_color('T', alpha)
    C_color = get_nuc_color('C', alpha)
    G_color = get_nuc_color('G', alpha)
    INDEL_color = get_nuc_color('N', alpha)

    if custom_colors is not None:
        hex_alpha = '66'  # this is equivalent to 40% in hexadecimal
        if 'A' in custom_colors:
            A_color = custom_colors['A'] + hex_alpha
        if 'T' in custom_colors:
            T_color = custom_colors['T'] + hex_alpha
        if 'C' in custom_colors:
            C_color = custom_colors['C'] + hex_alpha
        if 'G' in custom_colors:
            G_color = custom_colors['G'] + hex_alpha
        if 'N' in custom_colors:
            INDEL_color = custom_colors['N'] + hex_alpha

    dna_to_numbers = {'-': 0, 'A': 1, 'T': 2, 'C': 3, 'G': 4, 'N': 5}
    def seq_to_numbers(seq): return [dna_to_numbers[x] for x in seq]

    cmap = colors_mpl.ListedColormap(
        [INDEL_color, A_color, T_color, C_color, G_color, INDEL_color])

    # ref_seq_around_cut=reference_seq[max(0,cut_point-plot_nuc_len/2+1):min(len(reference_seq),cut_point+plot_nuc_len/2+1)]

#    print('per element anoot kws: ' + per_element_annot_kws)
    if len(per_element_annot_kws) > 1:
        per_element_annot_kws = np.vstack(per_element_annot_kws[::-1])
    else:
        per_element_annot_kws = np.array(per_element_annot_kws)
    ref_seq_hm = np.expand_dims(seq_to_numbers(reference_seq), 1).T
    ref_seq_annot_hm = np.expand_dims(list(reference_seq), 1).T

    annot = annot[::-1]
    X = X[::-1]

    N_ROWS = len(X)
    N_COLUMNS = plot_nuc_len

    if N_ROWS < 1:
        fig, ax = plt.subplots()
        fig.text(0.5, 0.5, 'No Alleles', horizontalalignment='center',
                 verticalalignment='center', transform=ax.transAxes)
        ax.set_clip_on(False)

        fig.savefig(fig_filename_root+'.pdf', bbox_inches='tight')
        if SAVE_ALSO_PNG:
            fig.savefig(fig_filename_root+'.png', bbox_inches='tight')
        plt.close(fig)
        return

    sgRNA_rows = []
    num_sgRNA_rows = 0

    if sgRNA_intervals and len(sgRNA_intervals) > 0:
        sgRNA_rows = get_rows_for_sgRNA_annotation(
            sgRNA_intervals, plot_nuc_len)
        num_sgRNA_rows = max(sgRNA_rows) + 1
        fig = plt.figure(
            figsize=(plot_nuc_len*0.3, (N_ROWS+1 + num_sgRNA_rows)*0.6))
        gs1 = gridspec.GridSpec(N_ROWS+2, N_COLUMNS)
        gs2 = gridspec.GridSpec(N_ROWS+2, N_COLUMNS)
        # ax_hm_ref heatmap for the reference
        ax_hm_ref = plt.subplot(gs1[0:1, :])
        ax_hm = plt.subplot(gs2[2:, :])
    else:
        fig = plt.figure(figsize=(plot_nuc_len*0.3, (N_ROWS+1)*0.6))
        gs1 = gridspec.GridSpec(N_ROWS+1, N_COLUMNS)
        gs2 = gridspec.GridSpec(N_ROWS+1, N_COLUMNS)
        # ax_hm_ref heatmap for the reference
        ax_hm_ref = plt.subplot(gs1[0, :])
        ax_hm = plt.subplot(gs2[1:, :])

    custom_heatmap(ref_seq_hm, annot=ref_seq_annot_hm, annot_kws={
                   'size': 16}, cmap=cmap, fmt='s', ax=ax_hm_ref, vmin=0, vmax=5, square=True)
    custom_heatmap(X, annot=np.array(annot), annot_kws={
                   'size': 16}, cmap=cmap, fmt='s', ax=ax_hm, vmin=0, vmax=5, square=True, per_element_annot_kws=per_element_annot_kws)

    ax_hm.yaxis.tick_right()
    ax_hm.yaxis.set_ticklabels(y_labels[::-1], rotation=True, va='center')
    ax_hm.xaxis.set_ticks([])

    if sgRNA_intervals and len(sgRNA_intervals) > 0:
        this_sgRNA_y_start = -1*num_sgRNA_rows
        this_sgRNA_y_height = num_sgRNA_rows - 0.3
        add_sgRNA_to_ax(ax_hm_ref, sgRNA_intervals, sgRNA_y_start=this_sgRNA_y_start, sgRNA_y_height=this_sgRNA_y_height, amp_len=plot_nuc_len,
                        font_size='small', clip_on=False, sgRNA_names=sgRNA_names, sgRNA_mismatches=sgRNA_mismatches, x_offset=0, label_at_zero=True, sgRNA_rows=sgRNA_rows)

# todo -- add sgRNAs below reference plot
#    if sgRNA_intervals:
#        ax_hm_anno=plt.subplot(gs3[2, :])
#        sgRNA_y_start = 0.3
# sgRNA_y_height = 0.1
#        sgRNA_y_height = 10
#        min_sgRNA_x = None
#        for idx,sgRNA_int in enumerate(sgRNA_intervals):
#            ax_hm_anno.add_patch(
#                patches.Rectangle((2+sgRNA_int[0], sgRNA_y_start), 1+sgRNA_int[1]-sgRNA_int[0], sgRNA_y_height,facecolor=(0,0,0,0.15))
#                )
#            #set left-most sgrna start
#            if not min_sgRNA_x:
#                min_sgRNA_x = sgRNA_int[0]
#            if sgRNA_int[0] < min_sgRNA_x:
#                min_sgRNA_x = sgRNA_int[0]
#        ax_hm_anno.text(2+min_sgRNA_x,sgRNA_y_start + sgRNA_y_height/2,'sgRNA ',horizontalalignment='right',verticalalignment='center')

    # print lines

    # create boxes for ins
    for idx, lss in insertion_dict.items():
        for ls in lss:
            ax_hm.add_patch(patches.Rectangle(
                (ls[0], N_ROWS-idx-1), ls[1]-ls[0], 1, linewidth=3, edgecolor='r', fill=False))

    # cut point vertical line
    if plot_cut_point:
        ax_hm.vlines([cutting_point], *ax_hm.get_ylim(), linestyles='dashed')

    ax_hm_ref.yaxis.tick_right()
    ax_hm_ref.xaxis.set_ticks([])
    ax_hm_ref.yaxis.set_ticklabels(['Reference'], rotation=True, va='center')

    gs2.update(left=0, right=1, hspace=0.05, wspace=0,
               top=1*(((N_ROWS)*1.13))/(N_ROWS))
    gs1.update(left=0, right=1, hspace=0.05, wspace=0,)

    sns.set_context(rc={'axes.facecolor': 'white', 'lines.markeredgewidth': 1,
                    'mathtext.fontset': 'stix', 'text.usetex': True, 'text.latex.unicode': True})

    proxies = [matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='black',
                                       mec='none', marker=r'$\mathbf{{{}}}$'.format('bold'), ms=18),
               matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='none',
                                       mec='r', marker='s', ms=8, markeredgewidth=2.5),
               matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='none',
                                       mec='black', marker='_', ms=2,)]
    descriptions = ['Substitutions', 'Insertions', 'Deletions']

    if plot_cut_point:
        proxies.append(
            matplotlib.lines.Line2D([0], [1], linestyle='--', c='black', ms=6))
        descriptions.append('Predicted cleavage position')

    # ax_hm_ref.legend(proxies, descriptions, numpoints=1, markerscale=2, loc='center', bbox_to_anchor=(0.5, 4),ncol=1)
    lgd = ax_hm.legend(proxies, descriptions, numpoints=1, markerscale=2,
                       loc='upper center', bbox_to_anchor=(0.5, 0), ncol=1, fancybox=True, shadow=False)

    fig.savefig(fig_filename_root+'.pdf', bbox_inches='tight',
                bbox_extra_artists=(lgd,))
    if SAVE_ALSO_PNG:
        fig.savefig(fig_filename_root+'.png', bbox_inches='tight',
                    bbox_extra_artists=(lgd,))
    plt.close(fig)


def plot_alleles_table(reference_seq, df_alleles, fig_filename_root, custom_colors, MIN_FREQUENCY=0.5, MAX_N_ROWS=100, SAVE_ALSO_PNG=False, plot_cut_point=True, sgRNA_intervals=None, cutting_point = 26.5, sgRNA_names=None, sgRNA_mismatches=None, annotate_wildtype_allele='****', **kwargs):
    """
    plots an allele table for a dataframe with allele frequencies
    input:
    reference_seq: the reference amplicon sequence to plot
    df_alleles: merged dataframe (should include columns "#Reads','%Reads')
    fig_filename: figure filename to plot (not including '.pdf' or '.png')
    MIN_FREQUENCY: sum of alleles % must add to this to be plotted
    MAX_N_ROWS: max rows to plot
    SAVE_ALSO_PNG: whether to write png file as well
    plot_cut_point: if false, won't draw 'predicted cleavage' line
    sgRNA_intervals: locations where sgRNA is located
    sgRNA_mismatches: array (for each sgRNA_interval) of locations in sgRNA where there are mismatches
    sgRNA_names: array (for each sgRNA_interval) of names of sgRNAs (otherwise empty)
    custom_colors: dict of colors to plot (e.g. colors['A'] = (1,0,0,0.4) # red,blue,green,alpha )
    annotate_wildtype_allele: string to add to the end of the wildtype allele (e.g. ** or '')
    """
    X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = prep_alleles_table(
        df_alleles, reference_seq, MAX_N_ROWS, MIN_FREQUENCY)
    if annotate_wildtype_allele != '':
        for ix, is_ref in enumerate(is_reference):
            if is_ref:
                y_labels[ix] += annotate_wildtype_allele
    plot_alleles_heatmap(reference_seq, fig_filename_root, X, annot, y_labels, insertion_dict, per_element_annot_kws,
                         custom_colors, SAVE_ALSO_PNG, plot_cut_point, sgRNA_intervals, sgRNA_names, sgRNA_mismatches,cutting_point)
