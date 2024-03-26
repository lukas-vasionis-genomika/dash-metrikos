import os
from pprint import pprint
import seaborn as sns
from utils_local import st_objects as sto
from plotly import express as px
import streamlit as st
import matplotlib.pyplot as plt

fq_path = "inputs/fastq"
fq_files = [x for x in os.listdir("inputs/fastq") if x.startswith('control')]
fq_obj = [sto.MyFastQ(fastq_path=fq_path, fastq_name=fq, format_extension='fastq') for fq in fq_files
          ]

fq = fq_obj[0]

f5_path = "D:/duomenys is labes kompo/03_21_81bp_PCR/AAA444/20240321_1208_MN39294_ASO266_f6112776"

fq = fq.assign_fast5(fast5_path=f5_path)
fq=fq.get_reads_length_stats(rm_outliers=True, override_self_scores=True)
fq=fq.distribution_qscore_n_nt()


"""
graphs
"""

# Score by nt index BOX PLOT

# fig=px.box(data_frame=fq.scores_over_nt_index,x='nt_index_bin',y="nt_score", points=False)
# fig.update_layout(xaxis_type='category')
# fig.show()