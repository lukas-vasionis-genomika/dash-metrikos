import pandas as pd
import numpy as np
from matplotlib import gridspec

import matplotlib.pyplot as plt
def forming_pairs_of_sub_plots(bc_fq_files,gridspec):
    gs=gridspec

    ax_pairs = {}
    fq_f_names = [x.fastq_name for x in bc_fq_files]
    for fq_f_name in fq_f_names:
        ax_pairs[f'{fq_f_name}'] = {}

        fq_index = fq_f_names.index(fq_f_name)
        gs_pair_index = fq_index * 2

        ax_pairs[f'{fq_f_name}']['ax1'] = plt.subplot(gs[gs_pair_index])
        ax_pairs[f'{fq_f_name}']['ax2'] = plt.subplot(gs[gs_pair_index + 1])

    return ax_pairs


def get_msno_data(fq, score_treshold, seq_length_treshold=None):
    """
    Gets dataframe that is compatable input for msno graph. That is, rows=reads; cols nucleotides of respective read
    if score_treshold is defined, replaces nt's whose scores are below score_treschold (not inclusive) with None. Thus,
    the resulting dataframe helps to show location of low score nucleotides in each read.

    :param fq: MyFastq object that holds read data in fastq format
    :param score_treshold: False/None to display all nucleotides; int the minimum score of nucleotide quality that one
    wants to display.
    :param seq_length_treshold: filters the reads below this lenght. Has two purposes:
        1) Sets the amount of columns (nuceotides) to display in the msno graph
        2) Discards reads with
    :return: pd.DataFrame
    """

    data = sorted(fq.scores, key=len)
    # Removing reads below length treshold
    data = [x for x in data if len(x) <= seq_length_treshold]
    # Replacing nucleotides bellow score treshold with none (if score_treshold provided)
    if score_treshold:
        data = [[None if j < score_treshold else j for j in row] for row in data]

    data = [xi + [None] * (seq_length_treshold - len(xi)) for xi in data]

    df = pd.DataFrame(data, columns=[str(x) for x in np.arange(0, seq_length_treshold, 1)])
    return df


def get_labels(seq_length_treshold, step):
    """
    Creates list of labels and ticks for msno graph
    :param seq_length_treshold: (same as get_msno_data) display first n nucleotides of the reads (this number sets the amount of columns
    (nuceotides) to display in the msno graph
    :param step: interval betwean tich labels
    :return:
    """

    labels = [x for x in np.arange(0, seq_length_treshold, step)]
    # Adding None to create intervals between the ticks
    labels = [[x] + [None] * (step - 1) for x in labels]
    # Flattening the list of lists
    labels = [item for row in labels for item in row]
    return labels
