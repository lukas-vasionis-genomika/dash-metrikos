import os
from pprint import pprint

from Bio import SeqIO
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import numpy as np


class MyFastQ:
    def __init__(self, fastq_path, fastq_name, format_extension):

        self.path_fastq = fastq_path
        self.fastq_name = fastq_name
        self.experiment = fastq_name.split('.')[0]
        self.basecaller = fastq_name.split('.')[1]
        self.precessing_stage = fastq_name.split('.')[2]
        self.reads = SeqIO.parse(f"{self.path_fastq}/{self.fastq_name}", format=format_extension)
        self.scores = [record.letter_annotations["phred_quality"] for record in self.reads]
        self.scores_over_nt_index = None

        # read summary
        # see assign_fast5()
        self.fast5_summary = None
        self.read_count = None
        self.read_pass_count = None
        self.read_fail_count = None

        """
        read_stats
        """

        # Length
        # get_reads_length_stats()
        self.lengths = None
        self.length_median = None
        self.lengths_pass = None
        self.lengths_fail = None

        # Score
        self.score_avg_per_read = None
        self.score_mean_per_read = None

    def assign_fast5(self, fast5_path):
        summary_file = [x for x in os.listdir(fast5_path) if x.startswith("sequencing_summary")]

        if summary_file:
            summary_file = summary_file[0]

            df = pd.read_csv(f"{fast5_path}/{summary_file}", sep='\t')
            if df is not None:
                self.read_count = df.shape[0]

                df_pass_filtering = df.groupby('passes_filtering').size()
                self.read_pass_count = df_pass_filtering[True]
                self.read_fail_count = df_pass_filtering[False]

        else:

            print("In Fast5 directory: no sequencing_summary* file found")
            print("Here's what's inside the directory:")
            print(os.listdir(fast5_path))
            df = None

        self.fast5_summary = df
        return self

    def get_reads_length_stats(self, rm_outliers=True, override_self_scores=False):
        """
        Gets basic read length statistics. If selected, removes read outliers by length
        :param rm_outliers: if True, removes outliers by 1.5 IQR
        :param override_self_scores: If True the property of self.scores is overriden with list of filtered read scores (if rm_outliers=True)
        :return:
            self.lengths : list
            self.lengths_nedian int
                If rm_outliers:
            self.lengths_pass : list
            self.lengths_fail : list

                if override_scores:
            self.scores : list

        """
        lengths = [len(x) for x in self.scores]

        if rm_outliers:
            Q1 = np.percentile(lengths, 25, method='midpoint')
            Q3 = np.percentile(lengths, 75, method='midpoint')
            IQR = Q3 - Q1
            upper = Q3 + 1.5 * IQR
            lower = Q1 - 1.5 * IQR
            # print(upper)
            # print(lower)

            lengths_pass = [x for x in lengths if x > lower and x < upper]
            lengths_fail = [x for x in lengths if x < lower and x > upper]

            # lengths_fail_lower= np.where((lengths <= lower) | (lengths >= upper))
            self.lengths_pass = lengths_pass
            self.lengths_fail = lengths_fail

            if override_self_scores:
                self.scores = [x for x in self.scores if len(x) > lower and len(x) < upper]

        self.lengths = lengths
        self.length_median = np.median(lengths)

        return self

    def distribution_qscore_n_nt(self, set_bins_nt_index=True):
        """
        Transforms MyFastq.reads into long pd.DataFrame of nucleotide score per read, where columns=["read_id", 'nt', "nt_score"]
        :param fq:

        :return: pd.DataFrame(columns nt, nt_score)
        """

        data = self.scores
        data = [[(nt_index + 1, nt_score) for nt_index, nt_score in enumerate(read)] for read in data]

        data = [x for xs in data for x in xs]

        df_data = pd.DataFrame(data, columns=['nt', "nt_score"])
        df_data = df_data.sort_values(by='nt')

        if set_bins_nt_index:
            nt_index = np.array(df_data["nt"])

            conditions = [
                (1 <= nt_index) & (nt_index < 5),
                (5 <= nt_index) & (nt_index < 10),
                (10 <= nt_index) & (nt_index < 15),
                (15 <= nt_index) & (nt_index < 20),
                (20 <= nt_index) & (nt_index < 30),
                (30 <= nt_index) & (nt_index < 40),
                (40 <= nt_index) & (nt_index < 50),
                (50 <= nt_index) & (nt_index < 60),
                (60 <= nt_index) & (nt_index < 70),
                (70 <= nt_index) & (nt_index < 80),
                (80 <= nt_index) & (nt_index < 90),
                (90 <= nt_index) & (nt_index < 100),
                (100 <= nt_index) & (nt_index < 120),
                (120 <= nt_index) & (nt_index < 140),
                (140 <= nt_index) & (nt_index < 160),
                (160 <= nt_index) & (nt_index < 200),
                (200 <= nt_index) & (nt_index < 250),
                (250 <= nt_index) & (nt_index < 300),
                (300 <= nt_index)]
            choices = ['1-4', '5-9', '10-14', '15-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70-79', '80-89',
                       '90-99', '100-119', '120-139', '140-159', '160-199', '200-249', '250-299', '300+']

            nt_index = np.select(conditions, choices)
            df_data.loc[:, 'nt_index_bin'] = nt_index

            df_data['nt_index_bin'] = df_data['nt_index_bin'].astype(str)

        self.scores_over_nt_index = df_data
        return self

    def get_reads_score_stats(self):

        self.score_mean_per_read = [np.mean(x) for x in self.scores]
        self.score_median_per_read = [np.median(x) for x in self.scores]
        return self
