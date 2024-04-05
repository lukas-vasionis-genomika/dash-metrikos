import os

import streamlit as st
from utils_local import utils as lu
from utils_local import objs_batteries as obj_b
from utils_local import objs_singles as obj_s

path_graph_objects="outputs/graph_objects"
for i in os.listdir(path_graph_objects):
    print(i)

st.markdown("## Read counts")
fig=lu.BatteryFig(path=f"{path_graph_objects}/fig_read_count_hist_F5").load_fig()
st.plotly_chart(fig)

st.markdown("## Read lengths")
fig=lu.BatteryFig(path=f"{path_graph_objects}/fig_read_length_box").load_fig()
st.plotly_chart(fig)
fig=lu.BatteryFig(path=f"{path_graph_objects}/fig_read_length_hist").load_fig()
st.plotly_chart(fig)

st.markdown("## Read scores")
fig=lu.BatteryFig(path=f"{path_graph_objects}/fig_score_box").load_fig()
st.plotly_chart(fig)
fig=lu.BatteryFig(path=f"{path_graph_objects}/fig_fq_score_cum_percent").load_fig()
st.plotly_chart(fig)

# fig=lu.BatteryFig(path=f"{path_graph_objects}/fig_scores_over_nt_index").load_fig()
# st.plotly_chart(fig)

st.markdown("## Occurences")

fig = lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_occurences_primer_counts_stack').load_fig()
st.plotly_chart(fig)
# fig = lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_occurences_primer_counts_group').load_fig()
# st.plotly_chart(fig)
fig = lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_occurences_fragment_counts_stack').load_fig()
st.plotly_chart(fig)
# fig = lu.BatteryFig(fig=fig, path='outputs/graph_objects/fig_occurences_fragment_counts_group').load_fig()
# st.plotly_chart(fig)

"""
Lentel4
"""

path_root = "/home/lukasv/projects/sekoskaita/tyrimai/short_read_problem/"
battery_fq = obj_b.BatteryFQ().load_battery_fastq(path_root)
graph_target_seq_counts = obj_b.SubstringCounts(battery_fq)
graph_target_seq_counts = graph_target_seq_counts.get_data().data

df_primer=graph_target_seq_counts.loc[graph_target_seq_counts['seq_type']=="primer",:]
st.table(df_primer)
df_fragment=graph_target_seq_counts.loc[graph_target_seq_counts['seq_type']=="fragment",:]
st.table(df_fragment)