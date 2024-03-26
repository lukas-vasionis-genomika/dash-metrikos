import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
from Bio import SeqIO
import os
from plotly.subplots import make_subplots
from tqdm import tqdm
from utils import transforms_utils as tu
import numpy as np
from tslearn.utils import to_time_series_dataset

from utils import fastq_utils as fqu
from utils import ml_utils as mlu
from matplotlib import gridspec
import utils.my_msno.missingno as msno
import matplotlib.pyplot as plt
from utils.msno_read_graph import utils as msno_r
import my_streamlit as st
import pylab as pl
from utils import graph_utils as gu

def get_fq_score_cumulative_percentage(fastq_path):
    """
    Joins commands from jastq utils to extract data from reads and plot cummulative percentage of scores
    :param fastq_path: path to fastq files dir
    :return: shows fig_q_score_groups_by_basecaller_n_trim_levels express graph
    """

    fq_files = [f for f in os.listdir(fastq_path) if f != "fastqc"]

    # Extracting and trandsforming fastq data
    fq_files = [fqu.MyFastQ(fastq_path=fastq_path, fastq_name=fq, format_extension='fastq') for fq in fq_files]

    # Forming score counts into pd.DataFrame() from each fastq score data
    data = pd.concat([fq.get_fastq_score_counts() for fq in fq_files])
    # Sending concated pd.DataFrame to fig_q_score_groups_by_basecaller_n_trim_levels express plot
    fig = fqu.get_score_cumsum_counts_fig(data)
    fig.show()


def get_clusters_of_sequences(fq):
    """
    pipeline transforms location and scoro data of nucleotide into time series.
    This lets to aggregate reads into clusters and analyze them by clusters.
    However, the shiluete score failed to reach 0.5 score in many k values from 2 to 100 and to 900

    This means that either the data variates too much, or im picking wrong metrics (like dtw) or evaluation methods (like shiluette).

    :param fq: MyFastq class from local fastq_utils files
    :param assert_csv_read_nt_score:
    :return:
    """

    # Check if primary dataframe with scores is created and assigned to MyFastq.df_distribution_qscore_n_nt
    # If not: assign
    if fq.df_distribution_qscore_n_nt is None:
        fq.get_distribution_qscore_n_nt()

    # Additional transformations
    ml_data = {k: list(v) for k, v in fq.df_distribution_qscore_n_nt.groupby('read_id')['nt_score'].apply(list).items()}
    ml_data = [np.array(x) for x in ml_data.values()]
    formatted_dataset = to_time_series_dataset(ml_data)
    # Creating models clusters and labels
    mlu.get_number_of_cluster_statistics(formatted_dataset, fastq_name=fq.fastq_name)


def get_fig_q_score_groups_by_basecaller_n_trim_levels(fq_files, list_tresholds):
    """
    For each basecaller (bc) for each file in bc dir,
    Trims read ends on defined range

    :param bc_path:
    :param basecallers:
    :return:
    """

    # fq_files = []
    # for bc in basecallers:
    #     fastq_path = f"{bc_path}/{bc}"
    #
    #     bc_fq_files = [f for f in os.listdir(fastq_path) if f.endswith('.fastq')]
    #     bc_fq_files = [fqu.MyFastQ(fastq_path=fastq_path, fastq_name=fq, basecaller=bc, format_extension='fastq') for fq
    #                    in bc_fq_files]
    #     fq_files += bc_fq_files

    dfs_fig = []
    for fq in tqdm(fq_files):
        fq.extract_scores()

        # read_count=len(fq.scores)

        df = tu.get_data_min_scores(fq, list_tresholds=list_tresholds)
        df.astype({'trim_level': 'str'})
        dfs_fig.append(df)
    df = pd.concat(dfs_fig)
    # fig = px.histogram(df, y="len_score_sequences_read_trimmed", color="trim_i_min_scores",
    #                    barmode='group', nbins=20)
    # fig.show()

    df = df.sort_values(by=['trim_level', "file_name"], ascending=[False, True])
    fig = px.bar(df, x="adapter trimming", y="count", color='min_score',
                 facet_col="trim_level", facet_row='basecaller',
                 title=f'Counts of mininmum read q-score (after trimming ends at level "x")',
                 color_discrete_map=gu.my_colors,
                 category_orders={"min_score": ["None","Empty read", "..9", "10..19", "20..29", "30..39", "40..."],
                                  "trim_level": sorted(df["trim_level"].unique().tolist()),
                                  "file_name": sorted(df["file_name"].unique().tolist()),
                                  "adapter trimming": ["raw", "trimmed"],
                                  "basecaller": [fq.basecaller for fq in fq_files]
                                  },
    width=1000,
    height=400
    )
    fig.update_layout(
        font=dict(
            size=5,  # Set the font size here
        ),
        legend_traceorder="reversed"
    )

    return fig

def get_fig_q_score_groups_by_basecaller_n_trim_levels_lengths(fq_files, list_tresholds):
    """
    For each basecaller (bc) for each file in bc dir,
    Trims read ends on defined range

    :param bc_path:
    :param basecallers:
    :return:
    """

    # fq_files = []
    # for bc in basecallers:
    #     fastq_path = f"{bc_path}/{bc}"
    #
    #     bc_fq_files = [f for f in os.listdir(fastq_path) if f.endswith('.fastq')]
    #     bc_fq_files = [fqu.MyFastQ(fastq_path=fastq_path, fastq_name=fq, basecaller=bc, format_extension='fastq') for fq
    #                    in bc_fq_files]
    #     fq_files += bc_fq_files

    dfs_fig = []
    for fq in tqdm(fq_files):
        fq.extract_scores()

        # read_count=len(fq.scores)

        df = tu.get_data_min_scores_and_lengths(fq, list_tresholds=list_tresholds)
        df.astype({'trim_level': 'str'})
        dfs_fig.append(df)
    df = pd.concat(dfs_fig)
    # print(df.head().to_string())
    # df.to_csv("outputs/data/fig_q_score_groups_by_basecaller_n_trim_levels_lengths/01_29_RE_lig_dorado_hac_sup.raw.csv",
    #           index=False)

    fig = px.histogram(df, x="len_score_sequences_read_trimmed",
                       color="trim_i_min_scores",
                       color_discrete_map=gu.my_colors,
                       facet_col="trim_level", facet_row='basecaller',
                       category_orders={
                           "min_score": ["Empty read", "..9", "10..19", "20..29", "30..39", "40..."],
                           "trim_level": sorted(df["trim_level"].unique().tolist()),
                           "file_name": sorted(df["file_name"].unique().tolist()),
                           "adapter trimming": ["raw", "trimmed"],
                           "basecaller": [fq.basecaller for fq in fq_files]},
                       barmode="group"
                       )
    # fig.add_annotation(
    #
    #     text="An annotation referencing the axes",
    #     row=1,
    #     col=1
    # )
    fig.update_xaxes(maxallowed=200)
    fig.show()
    # df = df.sort_values(by=['trim_level', "file_name"], ascending=[False, True])
    # fig = px.bar(df, x="adapter trimming", y="count", color='min_score',
    #              facet_col="trim_level", facet_row='basecaller',
    #              title=f'Counts of mininmum read q-score (after trimming ends at level "x")',
    #              category_orders={"min_score": ["None","Empty read", "..9", "10..19", "20..29", "30..39", "40..."],
    #                               "trim_level": sorted(df["trim_level"].unique().tolist()),
    #                               "file_name": sorted(df["file_name"].unique().tolist()),
    #                               "adapter trimming": ["raw", "trimmed"],
    #                               "basecaller": basecallers
    #                               },
    #              # width=1000,
    #              # height=400
    #              )
    fig.update_layout(
        font=dict(
            size=10,  # Set the font size here
        ),
        legend_traceorder="reversed"
    )
    return fig


def get_msno_read_graph(bc_fq_files, score_treshold, seq_length_treshold):
    """
    Creates matplotlib graph based on customised mnso library (see utils.my_msno).
    First the fastq file reads are filtered by sequence lenght:seq_length_treshold
    Then each fasqfile gets a subplot:
        1) left subplot represents raw reads (nucletide scores are not modified)
        2) right subplot represents reads where all nucleotide scores above score_treshold are turned into None
    Thus the produced graph shows how which nucleotides are of bad quality

    :param bc_fq_files: list of MyFastq objects with fastq_path, fastq_name and basecaller
    :return: msno.matrix fig
    """

    fig_msno_reads = plt.figure(figsize=(100, 300))
    gs = gridspec.GridSpec(len(bc_fq_files), 2, figure=fig_msno_reads, width_ratios=(15, 15))
    gs.update(wspace=0.08, hspace=0.5)

    ax_pairs = msno_r.forming_pairs_of_sub_plots(bc_fq_files, gridspec=gs)

    for fq in tqdm(bc_fq_files):
        # Reading each fastq file and extracting nt scores of its reads

        fq.extract_scores()

        read_count = len(fq.scores)

        msno.matrix(msno_r.get_msno_data(fq, seq_length_treshold=seq_length_treshold, score_treshold=None),
                    ax_title=f"{fq.fastq_name} All scores",
                    sparkline=False,
                    ax=ax_pairs[fq.fastq_name]['ax1'],
                    labels=msno_r.get_labels(seq_length_treshold, 20),
                    row_labels=np.arange(0, read_count - 1, step=3000),
                    fontsize=12
                    )

        msno.matrix(msno_r.get_msno_data(fq, seq_length_treshold=seq_length_treshold, score_treshold=score_treshold),
                    ax_title="Scores above threshold",
                    sparkline=False,
                    ax=ax_pairs[fq.fastq_name]['ax2'],
                    labels=msno_r.get_labels(seq_length_treshold, 20),
                    row_labels=np.arange(0, read_count - 1, step=3000),
                    fontsize=12
                    )

    # fig_msno_reads.savefig(f"outputs/msno_read_length_L{seq_length_treshold}_S{score_treshold}.png", format="png")
    return fig_msno_reads