import gzip

from Bio import SeqIO
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import utils.transforms_utils as tu
import numpy as np

class MyFastQ:
    def __init__(self, fastq_path, fastq_name, format_extension,basecaller=None):

        self.path_fastq = fastq_path
        self.fastq_name = fastq_name
        self.format_extension = format_extension
        self.reads = SeqIO.parse(f"{self.path_fastq}/{self.fastq_name}", format=format_extension)
        self.cum_score_counts=None # pd.Dataframe with score counts, cumsum, cumpercent.
        self.scores=None
        self.df_distribution_qscore_n_nt = None
        self.basecaller=basecaller

    def normalize_lengths(self, min_read_length=200):
        """
        Adds 0 to each read shorter than read_lemgth to make it at least as long as min_read_length
        :param sc:
        :return: list of normalized scores
        """
        if self.scores==None:
            print("Scores not extracted: RUN fq.extract_scores()")
            exit()

        sc=self.scores
        sc = sorted(sc, key=len)
        # elongating the reads shorte than min_read_length to length of min_read_length
        sc = [x + [0] * (min_read_length - len(x)) for x in sc]
        # Trimming all reads to min_read_length
        sc = [x[0:min_read_length] for x in sc]

        return sc

    def filter_reads_by_length(self, length_range):
        length_range=[*range(length_range[0], length_range[1], 1)]
        self.reads = [record for record in self.reads if len(record.seq) in length_range]

    def get_distribution_qscore_n_nt(self, load_csv=None, reduce_read_ids=False):

        """
        Transforms MyFastq.reads into long format pd.DataFrame of nucleotide score per read, where columns=["read_id", 'nt', "nt_score"]

        Read or create distribution_qscore_n_nt pd.Dataframe
        :param load_csv: None, True (default csv directory), file path (read other file than in default path, overwrite (if overwrite)
        :param reduce_read_ids:
        :return: pd.Dataframe
        """
        # Check if primary dataframe with scores is created
        # if yes, read
        # Load if true or path is set to load_csv
        if type(load_csv) is str and not 'overwrite':
            df = pd.read_csv(f"{load_csv}")
        # If true: load default. If file not found - create
        elif load_csv is True:
            try:
                df = pd.read_csv(f"outputs/{self.fastq_name}.read_nt_scores.csv")
            except FileNotFoundError:
                print(f"File not found: outputs/{self.fastq_name}.read_nt_scores.csv")
                print(f"Creating file: outputs/{self.fastq_name}.read_nt_scores.csv" )
                df = tu.distribution_qscore_n_nt(self, out_format="csv")
        # if NO, write and assign to df variable
        else:
            df = tu.distribution_qscore_n_nt(self, out_format="csv")

        # Mark if transform long read_id's into numeric ids
        if reduce_read_ids:
            df['read_id'] = df.groupby('read_id').ngroup()

        self.df_distribution_qscore_n_nt = df

    def extract_scores(self):
        """
        from self.reads extracts scores as list of lists (reads)

        :return: list of lists
        """
        self.scores = [record.letter_annotations["phred_quality"] for record in self.reads]

    def get_fastq_score_counts(self):
        """
        Calculates cumsum of all scores in the fastq file (does not distinguish between reads)
        :return: pd.DataFrame of cumulative sum and cumulative percentage of score counts
        """
        # Getting the list of scores and flatting it as each read has its own list of scores

        scores = [item for row in self.scores for item in row]

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
        score_count.loc[:, 'file_name'] = self.fastq_name.replace(f".{self.format_extension}", "")
        self.cum_score_counts=score_count
        return self.cum_score_counts

    def get_hist_read_lenght(self, ref_len):
        l_seqs = [len(x.seq) for x in self.reads]

        fig =px.histogram(
            title=f'Histogram: Read length  (len_ref={str(ref_len)})         File: {self.fastq_name}',
            x=l_seqs,
            range_x=[0,max(l_seqs)],
            width=800,
            height=400
        )
        if type(ref_len)==int:
            ref_len=[str(ref_len)]
        elif type(ref_len)==list:
            ref_len=[str(x) for x in ref_len]
        else:
            print(f"The ref_len should be list or integer. Current ref_len type is {type(ref_len)}")

        for r_l in ref_len:
            fig.add_vline(x=int(r_l), line_width=3, line_dash="dash", annotation_text=f"{r_l}")
        fig.update_layout(
            xaxis_title="Read length",
            yaxis_title="Read count by length")
        fig.show()

    def get_hist_read_sqore(self, nucleotide_range: list):
        """

        Args:
            first_n_nucleotides: list of zero based indexes of nucleotide range one wants to see

        Returns: px.fig

        """
        # trimming to keep selected range of indexes
        scores = [row[nucleotide_range[0]:nucleotide_range[1]] for row in self.scores]

        # flattening
        scores = [item for row in scores for item in row]

        fig =px.histogram(
            title=f'Scores of nucletides in range: {nucleotide_range[0]+1}-{nucleotide_range[1]+1}',
            x=scores,
            width=800,
            height=400
        )
        fig.update_yaxes(range=[0,100_000])
        fig.update_xaxes(range=[0, 60])

        fig.update_layout(
            xaxis_title="Nucleotide phread scores",
            yaxis_title="Nt count by scores",

        )
        fig.show()

def get_3d_line_scores_by_nt_read(df):
    fig = px.line_3d(df, x='nt', y='read_id', z='nt_score', line_group='read_id',
                     # range_x=[0,200],
                     # range_y=[0,10000],
                     # range_z=[0,60]
                        )
    fig.update_yaxes(showticklabels=False)
    fig.update_layout(width=1000,
                      height=450,
                      font_size=11,
                      scene_aspectmode='manual',
                      # scene_aspectratio=dict(x=100, y=1000, z=50),# regulates scene diameters (default - cube)
                      # # # scene_camera_eye=dict(x=1.45, y=1.45, z=1),
                      # # template="none"
                      )
    fig.show()
def get_score_cumsum_counts_fig(cum_score_counts):
    """
    NOTE:
        Must be outside the MyFastq class as this function can use data of mutilple FastQs
    USAGE:
    data = pd.concat([fq.get_fastq_score_counts() for fq in fq_files])
    fig = fqu.get_fig(data)

    :param cum_score_counts:
    :return:fig_q_score_groups_by_basecaller_n_trim_levels express fig
    """

    fig = px.line(cum_score_counts, x="score", y="cum_score_count_percent", color="file_name", width=800, height=400)

    fig.add_hline(y=0.9,
                  annotation_text=f"90%",
                  annotation_position="top right")
    fig.add_hline(y=0.7,
                  annotation_text=f"70%",
                  annotation_position="top right")
    fig.add_hline(y=0.5,
                  annotation_text=f"50%",
                  annotation_position="top right")

    return fig

def get_freq_tbl_read_len(MyFastq):
    len_reads=[len(str(r.seq)) for r in MyFastq.reads]
    len_reads=pd.DataFrame(len_reads, columns=['read_length']).value_counts().reset_index()
    len_reads_sorted=len_reads.sort_values(by="read_length", ascending=False)
    print(len_reads_sorted.head(100).to_string())

