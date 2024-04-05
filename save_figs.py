import os
from pprint import pprint

import pandas as pd
import seaborn as sns
from utils_local import objs_singles as obj_s
from utils_local import objs_batteries as obj_b
from utils_local import utils as lu
from plotly import express as px
import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go

import plotly.io as pio  # keep this!

pio.templates.default = "plotly"  # keep this!

list_f5_paths = [
    "/home/lukasv/projects/sekoskaita/tyrimai/short_read_problem/12_18_PGR_RE/AAA444",
    "/home/lukasv/projects/sekoskaita/tyrimai/short_read_problem/01_29_RE_lig/AAA444",
    "/home/lukasv/projects/sekoskaita/tyrimai/short_read_problem/03_21_81bp_PCR/AAA444",
]
path_root = "/home/lukasv/projects/sekoskaita/tyrimai/short_read_problem/"
battery_fq = obj_b.BatteryFQ().load_battery_fastq(path_root)
"""
Graphs
"""

# # READ COUNT
# # From F5 summary
# fig=obj_b.BatteryF5([obj_s.MyFast5(x).get_summary() for x in list_f5_paths]).read_count_hist()
# fig=lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_read_count_hist_F5').save_fig()
# # From FQ
# fig= battery_fq.fig_read_count_hist()
# fig=lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_read_count_hist_FQ').save_fig()

# LENGTHS
# Hist
fig=obj_b.fig_read_length_hist(battery_fq,rm_outliers=500)
fig=lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_read_length_hist').save_fig()
exit()
# Box
fig=obj_b.fig_read_length_box(battery_fq,rm_outliers=500)
fig=lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_read_length_box').save_fig()

# SCORES
# Median score box plot
fig = obj_b.fig_score_box(battery_fq)
fig=lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_score_box').save_fig()
# Cumulative score count
fig=obj_b.fig_fq_score_cum_percent(battery_fq)
fig=lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_fq_score_cum_percent').save_fig()
# Scores over NT index
fig = obj_b.fig_scores_over_nt_index(battery_fq)
fig = lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_scores_over_nt_index').save_fig()


# Seq counts
graph_target_seq_counts = obj_b.SubstringCounts(battery_fq)
graph_target_seq_counts = graph_target_seq_counts.get_data()

fig = graph_target_seq_counts.get_fig_primer('stack')
fig = lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_occurences_primer_counts_stack').save_fig()
fig = graph_target_seq_counts.get_fig_primer('group')
fig = lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_occurences_primer_counts_group').save_fig()

fig = graph_target_seq_counts.get_fig_fragments('stack')
fig = lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_occurences_fragment_counts_stack').save_fig()
fig = graph_target_seq_counts.get_fig_fragments('group')
fig = lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_occurences_fragment_counts_group').save_fig()
exit()










_

"""
graphs
"""

# Score by nt index BOX PLOT

# fig=px.box(data_frame=fq.scores_over_nt_index,x='nt_index_bin',y="nt_score", points=False)
# fig.update_layout(xaxis_type='category')
# fig.show()
