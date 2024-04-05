import os
from pprint import pprint

import Bio
import pandas as pd
import seaborn as sns
from Bio import SeqIO

from utils_local import objs_singles as obj_s
from utils_local import objs_batteries as obj_b
from utils_local import utils as lu
from plotly import express as px
import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go

import plotly.io as pio  # keep this!

from utils_ciurlionis import transforms as cut

# path_twist_tb_reads = '/home/lukasv/projects/twist/himnas/12_01_Twist/AAA444/lvasionis/himnas_reads.txt'
path_twist_tb_reads = '/home/lukasv/projects/twist/ont_raw_trimmed.txt'

with open(path_twist_tb_reads, 'r') as f:
    f = f.readlines()
    f = [r.replace('\n', '') for r in f]
    f=[len(r) for r in f
       if len(r)<500
       ]
df = pd.DataFrame(f, columns=['reads'])
read_count=df.shape[0]
df = df.value_counts().reset_index().sort_values(by=['count'], ascending=[False])
print(df.shape)
# df = df[df['count']>1000]
print(df['count'].sum()/read_count)
print(df.to_string())
print(read_count)

fig = px.bar(df, y='count', x='reads'
             # labels={'x':'occurences of ',
             #                        'y':'count'}
             )
fig.show()
exit()

root = "/home/lukasv/projects/sekoskaita/tyrimai/short_read_problem"
experiments = [x for x in os.listdir(root) if x != 'other']

list_fq = []
for ex in experiments:
    fastq_path = f"{root}/{ex}/new_basecalls/BC/dorado_sup/"
    fastq_name = os.listdir(fastq_path)[0]
    fq = obj_s.MyFastQ(fastq_path, fastq_name, 'fastq')
    list_fq.append(fq)

battery_fq = obj_b.BatteryFQ(list_fq)

df = obj_b.SubstringCounts(battery_fq).get_data().data

df_primer = df[df['seq_type'] == 'primer']

df_fragment = df.loc[df['seq_type'] == 'fragment', [x for x in df.columns if x != "primer_type"]]

df_primer_pivot = df_primer.pivot_table(index=['sub_string', 'side', 'seq_type', 'primer_type', 'seq_transformation'],
                                        columns='experiment',
                                        values='count_occurences',
                                        aggfunc='sum').reset_index()

df_fragment_pivot = df_fragment.pivot_table(index=['sub_string', 'side', 'seq_type', 'seq_transformation'],
                                            columns='experiment',
                                            values='count_occurences',
                                            aggfunc='sum').reset_index()

# Fill NA values with 0 or another placeholder if necessary
# pivoted_df = pivoted_df.fillna(0)

print(df_primer_pivot.to_string())
print(df_fragment_pivot.to_string())
print(df.to_string())
# for ex in df['experiment'].unique():
#     df_ex = df.loc[df["experiment"] == ex, :]
#
#     print(df_ex.to_string())
