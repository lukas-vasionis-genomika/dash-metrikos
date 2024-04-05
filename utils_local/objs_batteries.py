import os

import Bio
import pandas as pd
from plotly import express as px
from utils_local import objs_singles as obj_s
import os
import numpy as np
import plotly.graph_objects as go


class BatteryFQ:
    def __init__(self, list_fq=None):
        self.list_fq = list_fq

    def load_battery_fastq(self, path_root):
        """
        Loads fastq files and froms a battery object from them
        :param path_root: location where the experiment directories are. Usualy: "/home/lukasv/projects/sekoskaita/tyrimai/short_read_problem/"
        :return: battery object
        """
        list_experiments = [x for x in os.listdir(path_root) if x != "other"]
        list_fq_names = [os.listdir(f"{path_root}/{experiment}/new_basecalls/BC/dorado_sup") for experiment in
                         list_experiments]
        list_fq_names = [item for row in list_fq_names for item in row]
        list_fq_paths = [f"{path_root}/{experiment}/new_basecalls/BC/dorado_sup" for experiment in list_experiments]

        list_fq = list(zip(list_fq_paths, list_fq_names))
        self.list_fq = [obj_s.MyFastQ(fastq_path=x[0], fastq_name=x[1], format_extension='fastq') for x in list_fq]
        return self

    def fig_read_count_hist(self):
        def form_data():
            all_rows = []

            for fq in self.list_fq:
                r_pass = fq.read_pass_count

                f5_rows = [[fq.experiment, "pass", r_pass]]
                all_rows += f5_rows

            df = pd.DataFrame(data=all_rows, columns=["experiment", "read_status_pre_BC", "read_count"])
            return df

        def form_fig(fig_data):
            fig = px.bar(data_frame=fig_data, x="experiment", y="read_count",
                         color="read_status_pre_BC",
                         color_discrete_map={"failed": "red", "pass": "slateblue"},
                         barmode='stack',
                         text_auto=True)
            return fig

        fig_data = form_data()
        fig = form_fig(fig_data)
        return fig


class BatteryF5:
    def __init__(self, list_f5):
        self.list_f5 = list_f5

    def read_count_hist(self):
        def form_data():
            all_rows = []

            for f5 in self.list_f5:
                experiment = f5.fast5_summary['experiment_id'].unique()[0]
                r_pass = f5.read_pass_count
                r_fail = f5.read_fail_count

                f5_rows = [[experiment, "pass", r_pass], [experiment, "failed", r_fail]]
                all_rows += f5_rows

            df = pd.DataFrame(data=all_rows, columns=["experiment", "read_status_pre_BC", "read_count"])
            return df

        def form_fig(fig_data):
            fig = px.bar(data_frame=fig_data, x="experiment", y="read_count",
                         color="read_status_pre_BC",
                         color_discrete_map={"failed": "red", "pass": "slateblue"},
                         barmode='stack',
                         text_auto=True)
            fig.update_layout(
                title_text="Read counts in experiments",
                xaxis_title="Experiment",
                yaxis_title="Read count",
            )
            return fig

        fig_data = form_data()
        fig = form_fig(fig_data)

        return fig


def fig_read_length_hist(battery_fq, rm_outliers):
    def get_data():
        list_fq = battery_fq.list_fq
        list_fq = [fq.get_reads_length_stats(rm_outliers=rm_outliers, override_self_scores=False) for fq in list_fq]
        all_data = []
        for fq in list_fq:
            rows = [[fq.experiment, read_length] for read_length in fq.lengths]
            all_data += rows
        df = pd.DataFrame(all_data, columns=['experiment', "read_length"])
        return df

    def get_fig(df):
        fig = px.histogram(df, x="read_length", color="experiment", facet_row="experiment", width=600, height=1000)
        fig.update_traces(opacity=0.50)
        fig.update_layout(
            title_text="Read length histogram over experiments",
            xaxis_title="Nucleotide index",
            yaxis_title="Nucleotide scores",
        )
        return fig

    df = get_data()
    fig = get_fig(df)

    return fig


def fig_read_length_box(battery_fq, rm_outliers):
    def get_data():
        list_fq = battery_fq.list_fq
        list_fq = [fq.get_reads_length_stats(rm_outliers=rm_outliers, override_self_scores=False) for fq in list_fq]
        all_data = []
        for fq in list_fq:
            rows = [[fq.experiment, read_length] for read_length in fq.lengths]
            all_data += rows
        df = pd.DataFrame(all_data, columns=['experiment', "read_length"])
        return df

    def get_fig(df):
        fig = px.box(df, y="read_length", x="experiment", points=False)
        fig.update_traces(opacity=0.50)
        fig.update_layout(
            title_text="Lengths of reads",
            xaxis_title="Experiment",
            yaxis_title="Read length",
        )
        return fig

    df = get_data()
    fig = get_fig(df)

    return fig


def fig_score_box(battery_fq):
    def get_data_stats():
        data_all = {}
        for fq in battery_fq.list_fq:
            scores = fq.scores
            fq_scores_flat = [x for read in scores for x in read]

            fq_stats = [
                min(fq_scores_flat),
                np.percentile(fq_scores_flat, 25, method='midpoint'),
                np.median(fq_scores_flat),
                np.percentile(fq_scores_flat, 75, method='midpoint'),
                max(fq_scores_flat)
            ]

            data_all[fq.experiment] = fq_stats

        return data_all

    def get_fig_stats(data):
        fig = go.Figure()

        for experiment, stats in data.items():
            fig.add_trace(go.Box(
                x=[experiment],  # Assign experiment name to the x-axis
                lowerfence=[stats[0]],
                q1=[stats[1]],
                median=[stats[2]],
                q3=[stats[3]],
                upperfence=[stats[4]],
                name=experiment  # Box plot name
            ))

        fig.update_layout(
            title_text="Score distribution for each experiment",
            xaxis_title="Experiment",
            yaxis_title="Score",
            boxmode='group'  # Groups box plots of the same x value next to each other
        )
        return fig

    data = get_data_stats()
    fig = get_fig_stats(data)
    return fig


def fig_fq_score_cum_percent(battery_fq):
    """
    Joins commands from jastq utils to extract data from reads and plot cummulative percentage of scores
    :param fastq_path: path to fastq files dir
    :return: shows fig_q_score_groups_by_basecaller_n_trim_levels express graph
    """

    def get_data():
        # Extracting and trandsforming fastq data
        fq_files = battery_fq.list_fq
        data_all = []
        for fq in fq_files:
            """
            Calculates cumsum of all scores in the fastq file (does not distinguish between reads)
            :return: pd.DataFrame of cumulative sum and cumulative percentage of score counts
            """
            # Getting the list of scores and flatting it as each read has its own list of scores

            scores = [item for row in fq.scores for item in row]

            # Turning into pandas df, getting the counts of all scores in it
            score_count = pd.Series(scores).value_counts().reset_index()
            # Assigning columns names
            score_count.columns = ['score', "count"]

            # Sorting it for cumulative sum
            score_count = score_count.sort_values('score').reset_index(drop=True)
            # Calculating cumulative sum
            score_count.loc[:, 'cum_score_count'] = score_count["count"].cumsum()
            # Counting the percentage of score occurances
            score_count.loc[:, 'cum_score_count_percent'] = score_count["cum_score_count"] / score_count["count"].sum()

            # Assigning the file name, removing extention from the file name
            score_count.loc[:, 'experiment'] = fq.experiment

            data_all.append(score_count)

        # Forming score counts into pd.DataFrame() from each fastq score data
        df = pd.concat(data_all)
        return df

    def get_fig(df):
        """
        NOTE:
            Must be outside the MyFastq class as this function can use data of mutilple FastQs
        USAGE:
        data = pd.concat([fq.get_fastq_score_counts() for fq in fq_files])
        fig = fqu.get_fig(data)

        :param cum_score_counts:
        :return:fig_q_score_groups_by_basecaller_n_trim_levels express fig
        """

        fig = px.line(df, x="score", y="cum_score_count_percent", color="experiment", width=800,
                      height=400)

        fig.add_hline(y=0.9,
                      annotation_text=f"90%",
                      annotation_position="top right")
        fig.add_hline(y=0.7,
                      annotation_text=f"70%",
                      annotation_position="top right")
        fig.add_hline(y=0.5,
                      annotation_text=f"50%",
                      annotation_position="top right")

        fig.update_layout(
            title_text="Cumulative score percent of all reads",
            xaxis_title="Score",
            yaxis_title="Cumulative percent",
        )
        return fig

    df = get_data()

    fig = get_fig(df)

    return fig


def fig_scores_over_nt_index(battery_fq):
    def get_data():
        list_fq = battery_fq.list_fq
        list_fq = [fq.distribution_qscore_n_nt() for fq in list_fq]
        all_data = [fq.scores_over_nt_index for fq in list_fq]
        df = pd.concat(all_data)
        return df

    def get_fig(df):
        fig = px.box(data_frame=df, x='nt_index_bin', y="nt_score", color='experiment', points=False)
        fig.update_layout(xaxis_type='category')

        fig.update_layout(
            title_text="Nucleotide score distribution over read nucleotides",
            xaxis_title="Nucleotide index",
            yaxis_title="Nucleotide scores",
        )
        return fig

    df = get_data()
    figure = get_fig(df)

    return figure


class SubstringCounts:
    def __init__(self, battery_fq, data=None, fig=None):
        self.battery_fq = battery_fq
        self.data = data
        self.fig = fig

    def get_data(self, substring_dict=None, group_by_transformation=False):
        """
        Extracts occurences of substrings from substring dict and puts into dataframe
        The returned dataframe is compatable with px.Bar() plot
        :param group_by_transformation:
        :param substring_dict:
        :return:
        """

        # Suplemental functions
        def supplement_dict_substring(target_sub_strings):
            for key in list(target_sub_strings.keys()):
                value = Bio.Seq.Seq(target_sub_strings[key])
                # Reverse the sequence
                target_sub_strings[f"{key}_R"] = str(value[::-1])

                # Complement the sequence
                target_sub_strings[f"{key}_C"] = str(value.complement())

                # Reverse complement the sequence
                target_sub_strings[f"{key}_RC"] = str(value.reverse_complement())
            return target_sub_strings

        def find_sub_string_occurences(strings, reads_sequences, experiment):
            # total_reads = len(reads_sequences)

            all_substring_rows = []
            for key, value in strings.items():
                contains_substring = list(map(lambda x: value in x, reads_sequences))
                occured = contains_substring.count(True)

                key_split = key.split("_")
                side = key_split[0]
                seq_type = key_split[1]
                primer_type = key_split[2] if seq_type == 'primer' else None
                transformation = key_split[-1] if key_split[-1] in ["R", "RC", "C"] else 'RAW'

                row = [experiment, key, side, seq_type, primer_type, transformation, occured]

                all_substring_rows.append(row)

            df = pd.DataFrame(all_substring_rows,
                              columns=["experiment", "sub_string", "side", "seq_type", "primer_type",
                                       "seq_transformation", "count_occurences"])


            return df

        """
        Data processing
        """

        list_fq = self.battery_fq.list_fq

        if substring_dict is None:
            print("Object SubstringCounts: Using default dictionary of fragment and primer sequences")
            sub_strings = {
                'top_primer_F': "CTACAACGCAGATTACAACCTCAGTG",
                'top_fragment': "CCCACAAACATGCCCAACAACATCACCAGCACCAAT",
                'top_primer_R': "GGTAACGCTGGCAAGGATAGGAA",
                'bot_primer_F': "AAGATGTTGCGTCTAATGTTGGAGTCAC",
                'bot_fragment': "GGGTGTTTGTACGGGTTGTTGTAGTGGTCGTGGTTA",
                'bot_primer_R': "CCATTGCGACCGTTCCTATCC"
            }
        else:
            sub_strings = substring_dict

        sub_strings = supplement_dict_substring(sub_strings)

        list_df = []
        for fq_obj in list_fq:
            seqs = fq_obj.extract_seqs().seqs

            df_fq = find_sub_string_occurences(sub_strings, seqs, fq_obj.experiment)
            list_df.append(df_fq)

        df_all = pd.concat(list_df)
        # df_all = df_all.loc[df_all['count_occurences'] != 0, :]
        df_all = df_all.sort_values(by=["experiment", "sub_string"],
                                    ascending=[False, False])
        self.data = df_all
        return self

    def get_fig_primer(self, barmode):
        df_primer = self.data.loc[self.data["seq_type"] == 'primer', :]

        if barmode == "group":
            fig = px.bar(df_primer, x="seq_transformation", y='count_occurences',
                         color='experiment', facet_col="primer_type", facet_row="side",
                         color_discrete_map={"03_21_81bp_PCR": '#00cc96', "01_29_RE_lig": "#EF553B",
                                             "12_18_PGR_RE": "#636efa"},
                         category_orders={'experiment': ["12_18_PGR_RE", "01_29_RE_lig", "03_21_81bp_PCR"],
                                          "seq_transformation": ['RAW', "R", "C", "RC"],
                                          "primer_type": ['F', "R"]
                                          },

                         height=700, barmode=barmode)

            fig.update_xaxes(title='')
            fig.update_yaxes(title='')

            fig.update_layout(title_text="Primer occurences in reads",
                              xaxis1=dict(
                                  title="Sequence Transformation <br> RAW, R:reverse, C:compliment, RC:reverse-complement"),
                              yaxis1=dict(title="Read count with single occurence"))

        elif barmode == "stack":
            fig = px.bar(df_primer, x="experiment", y='count_occurences',
                         color='seq_transformation', facet_col="primer_type", facet_row="side",
                         color_discrete_map={'RAW': "#000000", "R": "#869fda", "C": "#eea176", "RC": "#95bb88"},
                         category_orders={'experiment': ["12_18_PGR_RE", "01_29_RE_lig", "03_21_81bp_PCR"],
                                          "seq_transformation": ['RAW', "R", "C", "RC"],
                                          "primer_type": ['F', "R"]
                                          },
                         height=700, barmode=barmode)

            fig.update_xaxes(title='')
            fig.update_yaxes(title='')

            fig.update_layout(title_text="Primer occurences in reads",
                              xaxis1=dict(title="Experiment"),
                              yaxis1=dict(title="Read count with single occurence"))

        return fig

    def get_fig_fragments(self, barmode):

        df_fragment = self.data.loc[self.data["seq_type"] == 'fragment', :]
        if barmode == "group":
            fig = px.bar(df_fragment, x="seq_transformation", y='count_occurences',
                         color='experiment', facet_row="side",
                         color_discrete_map={"03_21_81bp_PCR": '#00cc96', "01_29_RE_lig": "#EF553B",
                                             "12_18_PGR_RE": "#636efa"},
                         category_orders={'experiment': ["12_18_PGR_RE", "01_29_RE_lig", "03_21_81bp_PCR"],
                                          "seq_transformation": ['RAW', "R", "C", "RC"]
                                          },
                         height=600, barmode=barmode)

            fig.update_xaxes(title='')
            fig.update_yaxes(title='')

            fig.update_layout(
                title_text="Fragment occurences in reads",
                xaxis1=dict(
                    title="Sequence Transformation <br> (RAW, R:reverse, C:compliment, RC:reverse-complement"),
                yaxis1=dict(title="Read count with single occurence"))

        elif barmode == "stack":
            fig = px.bar(df_fragment, x="experiment", y='count_occurences',
                         color='seq_transformation', facet_row="side",
                         color_discrete_map={'RAW': "#000000", "R": "#869fda", "C": "#eea176", "RC": "#95bb88"},
                         category_orders={'experiment': ["12_18_PGR_RE", "01_29_RE_lig", "03_21_81bp_PCR"],
                                          "seq_transformation": ['RAW', "R", "C", "RC"]
                                          },
                         height=600, barmode=barmode)

            fig.update_xaxes(title='')
            fig.update_yaxes(title='')

            fig.update_layout(
                title_text="Fragment occurences in reads",
                xaxis1=dict(title="Experiment"),
                yaxis1=dict(title="Read count with single occurence"))
        return fig
