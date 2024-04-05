from pprint import pprint

import pandas
import Bio
import pandas as pd

import utils_local.objs_singles as obj_s
import utils_local.objs_batteries as obj_b

fq_path = "/home/lukasv/projects/twist/03_31_ciurlionis"
fq_name = "03_31_ciurlionis.dorado_sup.raw.fastq"

fq = obj_s.MyFastQ(fq_path, fq_name, format_extension="fastq")

battery_fq = obj_b.BatteryFQ(list_fq=[fq])
s_counts = obj_b.SubstringCounts(battery_fq).get_data(substring_dict={"top_primer_5": "GTCAGATGTGTATAAGAGACAG",
                                                                      "bottom_primer_5 ": "CGTGTGCTCTTCCGATCT"},
                                                      group_by_transformation=True)

fig = battery_fq.fig_read_count_hist()
fig.show()
fig = obj_b.fig_read_length_box(battery_fq, rm_outliers=False)
fig.show()
fig = obj_b.fig_read_length_hist(battery_fq, rm_outliers=False)
fig.show()
print(s_counts.data)
exit()
# Išgaunu sekų id ir poras
seqs = [[r.id, str(r.seq)] for r in fq.reads]


# Suskaičiuoju kiek primerių randama sekoje
def find_substring_indices(s, sub):
    def find_indices(start, indices):
        index = s.find(sub, start)
        if index == -1:
            return [len(indices), indices]
        return find_indices(index + 1, indices + [index])

    return find_indices(0, [])


primer_5 = "GTCAGATGTGTATAAGAGACAG"
seqs = [r + find_substring_indices(r[1], primer_5) for r in seqs]

primer_5_RC = str(Bio.Seq.Seq(primer_5).reverse_complement())
seqs = [r + find_substring_indices(r[1], primer_5_RC) for r in seqs]

primer_3 = "CGTGTGCTCTTCCGATCT"
seqs = [r + find_substring_indices(r[1], primer_3) for r in seqs]

primer_3_RC = str(Bio.Seq.Seq(primer_3).reverse_complement())
seqs = [r + find_substring_indices(r[1], primer_3_RC) for r in seqs]

"""
Filter only those that have occurence
"""

P5_occurences = [x for x in seqs if x[2] != 0]
P5_RC_occurences = [x for x in seqs if x[4] != 0]
P3_occurences = [x for x in seqs if x[6] != 0]
P3_RC_occurences = [x for x in seqs if x[8] != 0]

print("P5 occurences", len([x[2] for x in P5_occurences]))
print("P5_RC occurences", len([x[4] for x in P5_RC_occurences]))
print("P3 occurences", len([x[6] for x in P3_occurences]))
print("P3_RC occurences", len([x[8] for x in P3_RC_occurences]))

"""
Getting examples
"""
# print("primer_5\t",primer_5)
# print(P5_occurences[0])
#
# print("primer_5_RC\t", primer_5_RC)
# print(P5_RC_occurences[0])
#
# print("primer_3\t", primer_3)
# print(P3_occurences[0])
#
# print("primer_3_RC\t", primer_3_RC)
# print(P3_RC_occurences[0])

"""
Toliau dirbsiu tik su readais kur aptiktas P3, nes jo daugiausiai rasta
"""

"""
Patikrinu, kur dažniausiai randami P3 pradmenys pradžioj, vidury ar gale read'o
"""


def transform_list(primary_list):
    # Unpack the first element and the sublist from the primary list
    first_element, sublist = primary_list[0], primary_list[1]
    # Use list comprehension to create the new list of lists
    transformed_list = [[first_element, elem] for elem in sublist]
    return transformed_list


reads_P3 = [transform_list([x[6], x[7]]) for x in P3_occurences]
reads_P3 = [x for row in reads_P3 for x in row]
df_r_P3 = pd.DataFrame(reads_P3, columns=['count_occurences', "occurence_position"])

import plotly.express as px


def show_P3_occurence_position_and_count(df_r_P3):
    fig = px.box(df_r_P3, x='count_occurences', y="occurence_position")
    fig.show()


# show_P3_occurence_position_and_count(df_r_P3)


"""
Daugelyje read'ų P3 primeriai  randasi viduryje arba gale. Vadinasi karpysiu nuo primerio į galą
STRATEGIJA:

Split reads and P3

Atlikti ilgių analizę kad nustatyti L

Kiekvienam read'ui imsiu paskutinius L nukleotidų

L=161, bet reikia pasirinkti ir trumpesnius vadovaujantis ilgiu analize
"""


# Ilgigių analizė sukarpius read'us
def length_analysis_after_split(P3_occurences):
    # spliting the reads
    reads_P3_cut = [r[1].split(primer_3) for r in P3_occurences]
    reads_P3_cut = [x for row in reads_P3_cut for x in row]

    # Calculating their lengths
    reads_P3_cut_len = [len(x) for x in reads_P3_cut]

    fig = px.histogram(x=reads_P3_cut_len)
    fig.show()


# Filtravimas pagal ilgius. Nustatytas L = 100
def filter_by_lenght(read_P3_cut, min_read_length):
    def filter_cut_reads(cut_reads):
        cut_reads = [x for x in cut_reads if len(x) >= min_read_length]
        return cut_reads

    read_P3_cut_filtered_len = [[r[0], filter_cut_reads(r[1])] for r in read_P3_cut]
    read_P3_cut_filtered_len = [r for r in read_P3_cut_filtered_len if r[1] != []]
    return read_P3_cut_filtered_len


# Forming the dataset
cut_read_set_to_filter = [[r[0], r[1].split(primer_3)] for r in P3_occurences]

# Filtering reads
read_P3_cut_filtered_len = filter_by_lenght(cut_read_set_to_filter, min_read_length=100)

seq_P3_cut_filtered_len = [x[1] for x in read_P3_cut_filtered_len]
seq_P3_cut_filtered_len = [x for row in seq_P3_cut_filtered_len for x in row]

# Read trimming - removing last L nucleotides from read
seq_P3_cut_filtered_len = [x[-100:] for x in seq_P3_cut_filtered_len]
print(len(seq_P3_cut_filtered_len))
"""
writing to a text
"""
# with open(f'outputs/{fq_name}_P3_minus_100.txt', mode='wt', encoding='utf-8') as myfile:
#     myfile.write('\n'.join(seq_P3_cut_filtered_len))
#     myfile.write('\n')


"""
Decoder:
root@8109721c6426:/# ls data/
twist.dorado_sup.raw.fastq_P3_minus_100.txt
root@8109721c6426:/# /usr/twist/bin/twist_import /data /data/imported
Found 1 read pools
Imported 25124 reads from /data/twist.dorado_sup.raw.fastq_P3_minus_100.txt
All done
root@8109721c6426:/# /usr/twist/bin/twist_cluster /data/imported /data/cluster
Building inner codec dom2
Finished building inner codec dom2
Found 1 read pools
Loading reads from twist.dorado_sup.raw.fastq_P3_minus_100.txt
Decoding read indexes
[####                                              ] 8%
Discarded 23045 reads
Clustering reads
[##################################################] 100%
Created 2058 clusters from 25124 reads
Calculating clustering statistics for read id twist.dorado_sup.raw.fastq_P3_minus_100.txt
Average cluster size = 1
0 clusters with source ids over 2058 total clusters
0 full matching clusters over 0 clusters with source ids
Finished writing clusters for twist.dorado_sup.raw.fastq_P3_minus_100.txt
All done
root@8109721c6426:/# /usr/twist/bin/twist_cluster_merge /data/imported /data/cluster /data/merged
Found 1 read pools
Loading reads from twist.dorado_sup.raw.fastq_P3_minus_100.txt
Aligning clusters from twist.dorado_sup.raw.fastq_P3_minus_100.txt
Generated 2058 aligned reads
All done
root@8109721c6426:/# /usr/twist/bin/twist_decode /data/merged /data/decoded
Using filesystem layer
Using greedy only method
Building codec with outer RS12 and inner dom2
Finished building codec
Found 1 reads
No clustering provided for twist.dorado_sup.raw.fastq_P3_minus_100.txt - creating single read clusters
Using 12 hardware concurrency
Running cluster inner decoding
[   #  #                #   #   #               #] 100%
Running frame outer decoding
[                                                  ] 0%terminate called after throwing an instance of 'std::invalid_argument'
  what():  Missing frame index 0 containing metadata
Aborted (core dumped)

"""
