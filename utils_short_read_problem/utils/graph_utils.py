import json
from pprint import pprint

import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
from Bio import SeqIO
import os
from plotly.subplots import make_subplots
from utils import fastq_utils as fqu
from utils import transforms_utils as tu
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import pickle
import io
from PIL import Image

my_colors = {"Empty read": "#000000",
             "..9": "#e29578",
             "10..19": "#FB5607",
             "20..29": "#8338EC",
             "30..39": "#3A86FF",
             "40...": "#3A86FF"}


class ONT_fig:
    def __init__(self,figure=None, f_path=None, data=None, d_path=None):
        self.data = data
        self.figure = figure
        self.f_path = f_path
        self.d_path = d_path

    def save_fig(self):
        pickle.dump(obj=self.figure,
                    file=open(self.f_path, "wb"))

    def load_fig(self):
        fig = pickle.load(
            file=open(self.f_path, "rb"))

        return fig

    def save_data(self, path):
        pass

    def load_data(self, path):
        pass

    def transform_data(self):
        pass
def save_fig_as_img(figure, path, format):

    figure.savefig(path, format=format)

def get_3d_line_scores_by_nt_read(fq):
    """
    REJECETED: too much data, the graph is too giberish. Use only for small fastq files
    :param fq: MyFastq class object
    :return:
    """
    if fq.df_distribution_qscore_n_nt is None:
        fq.get_distribution_qscore_n_nt()

    fig = px.line_3d(fq.df_distribution_qscore_n_nt, x='nt', y='read_id', z='nt_score', line_group='read_id',
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


def get_qc_score_graph(fq):
    fq.extract_scores()
    raw_scores = fq.scores
    nt_scores = [item for sublist in raw_scores for item in sublist]

    # Get indexes of respective elements in inner lists of l_x
    nt_index = [index + 1 for sublist in raw_scores for index, _ in enumerate(sublist)]

    g = sns.boxplot(y=nt_scores, x=nt_index)
    g.tick_params(axis='x', rotation=90)
    g.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=50))
    g.axhline(y=10)

    plt.show()

import numpy as np
def avg_score(fq_list):
    list_df=[]
    for fq in fq_list:
        fq.extract_scores()
        sc_avg = [np.mean(x) for x in fq.scores]

        df=pd.DataFrame({"score_avg":sc_avg})
        df.loc[:,"file_name"]=fq.fastq_name
        list_df.append(df)
    df_all=pd.concat(list_df)

    fig=px.histogram(df_all,x="score_avg", color="file_name", nbins=20, barmode="group")
    fig.show()
    return fig


