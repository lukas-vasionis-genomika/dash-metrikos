from Bio.Seq import Seq
import pandas as pd
from tqdm import tqdm
import plotly.express as px
import utils_ciurlionis.transforms as cut
from utils_local import objs_singles as obj_s
# seqs_file = "/home/lukasv/projects/twist/03_31_ciurlionis/03_31_ciurlionis.dorado_sup.raw.txt"

seqs_file = "/home/lukasv/projects/sekoskaita/tyrimai/short_read_problem/01_29_RE_lig/new_basecalls/BC/dorado_sup/01_29_RE_lig.dorado_sup.pass.raw.fastq"
experiment = seqs_file.split('/')[-1].split('.')[0]

substrings = {
    "target_RRS_1_RAW": "GGATG[A-Z][A-Z]",
    "target_RRS_1_C": "CCTAC[A-Z][A-Z]",
    "target_RRS_1_R": "[A-Z][A-Z]GTAGG",
    "target_RRS_1_RC": "[A-Z][A-Z]CATCC",
    # "target_RRS_2_RAW": "[A-Z][A-Z]CATCC",
    # "target_RRS_2_C ": "[A-Z][A-Z]GTAGG",
    # "target_RRS_2_R": "CCTAC[A-Z][A-Z]",
    # "target_RRS_2_RC": "GGATG[A-Z][A-Z]",
}

# with open(seqs_file, 'r') as seqs:
#     seqs = seqs.readlines()
#     seqs = [x.replace("\n", "") for x in seqs]
fastq_path="/home/lukasv/projects/sekoskaita/tyrimai/short_read_problem/12_18_PGR_RE/new_basecalls/BC/dorado_sup"
fastq_name="12_18_PGR_RE.dorado_sup.pass.raw.fastq"

fq=obj_s.MyFastQ(fastq_path=fastq_path, fastq_name=fastq_name, format_extension='fastq')
seqs=[str(x) for x in fq.extract_seqs().seqs]

s_counts = cut.get_count_occurences(seqs, substring_dict=substrings, supplement_dict_substring=False, is_regex=True)
s_counts.to_csv(f"outputs/ciurlionis/{fq.experiment}_counts_occurences_RRS.csv")

