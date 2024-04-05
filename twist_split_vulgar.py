from pprint import pprint

import pandas
import Bio
import pandas as pd

import utils_local.objs_singles as obj_s
import utils_local.objs_batteries as obj_b

fq_path = "/home/lukasv/projects/twist/himnas/12_01_Twist/AAA444/lvasionis"
fq_name = "twist.dorado_sup.raw.fastq"

fq = obj_s.MyFastQ(fq_path, fq_name, format_extension="fastq")

seqs = [str(r.seq) for r in fq.reads]
def split_string_at_len(string, lenght):

    res = [string[y - lenght:y] for y in range(lenght, len(string) + lenght, lenght)]
    return res

seqs = [split_string_at_len(x,100) for x in seqs]
seqs = [x for row in seqs for x in row]
seqs = [x for x in seqs if len(x)==100]

print(len(seqs))
"""
writing to a text
"""
with open(f'outputs/{fq_name}_l_100.txt', mode='wt', encoding='utf-8') as myfile:
    myfile.write('\n'.join(seqs))
    myfile.write('\n')

