from Bio.Seq import Seq
import pandas as pd
from tqdm import tqdm
import regex as re


def supplement_dict_substring(target_sub_strings):
    for key in list(target_sub_strings.keys()):
        value = Seq(target_sub_strings[key])
        # Reverse the sequence
        target_sub_strings[f"{key}_R"] = str(value[::-1])

        # Complement the sequence
        target_sub_strings[f"{key}_C"] = str(value.complement())

        # Reverse complement the sequence
        target_sub_strings[f"{key}_RC"] = str(value.reverse_complement())
    return target_sub_strings


def get_count_occurences(seqs, substring_dict=None, is_regex=False, supplement_dict_substring=True):
    """
    Extracts occurences of substrings from substring dict and puts into dataframe
    The returned dataframe is compatable with px.Bar() plot
    :param group_by_transformation:
    :param substring_dict:
    :return:
    """

    # Suplemental functions

    def find_sub_string_occurences(strings, reads_sequences, experiment):
        total_reads = len(reads_sequences)
        all_substring_rows = []
        for key, value in tqdm(strings.items()):
            if is_regex:
                contains_substring = list(map(lambda x: bool(re.search(value, x)), reads_sequences))
            else:
                contains_substring = list(map(lambda x: value in x, reads_sequences))
            occured = contains_substring.count(True)

            key_split = key.split("_")
            side = key_split[0]
            seq_type = key_split[1]
            primer_type = key_split[2] if seq_type == 'primer' else None
            transformation = key_split[-1] if key_split[-1] in ["R", "RC", "C"] else 'RAW'

            row = [experiment, key, side, seq_type, primer_type, transformation, occured]

            all_substring_rows.append(row)

        df = pd.DataFrame(all_substring_rows,
                          columns=["experiment", "sub_string", "side", "seq_type", "primer_type",
                                   "seq_transformation", "count_occurences"])
        df.loc[:, "total_reads"] = total_reads
        df.loc[:, "percentage_found"] = df["count_occurences"] / total_reads
        return df

    sub_strings = substring_dict
    if supplement_dict_substring:
        sub_strings = supplement_dict_substring(sub_strings)

    list_df = []

    df_fq = find_sub_string_occurences(sub_strings, seqs, "ciurlionis")
    list_df.append(df_fq)

    df_all = pd.concat(list_df)
    # df_all = df_all.loc[df_all['count_occurences'] != 0, :]
    df_all = df_all.sort_values(by=["experiment", "sub_string"],
                                ascending=[False, False])

    return df_all


def get_count_occurences_agg(seqs, substring_dict, is_regex=False):
    total_reads = len(seqs)
    data_contains = {}
    for key, value in tqdm(substring_dict.items()):
        if is_regex:
            contains_substring = list(map(lambda x: bool(re.search(value, x)), seqs))
        else:
            contains_substring = list(map(lambda x: value in x, seqs))
        data_contains[key] = contains_substring
    df = pd.DataFrame(data_contains)
    df = df.groupby(list(substring_dict.keys())).size().reset_index(name='count')

    df['sum_any'] = df[[x for x in df.columns if x != 'count']].sum(axis=1)

    df['total_reads'] = total_reads
    df['found_percent'] = df['count'] / total_reads
    #
    # for col in df.columns:
    #     df[col] = df[col].str.replace(False, col)

    return df
