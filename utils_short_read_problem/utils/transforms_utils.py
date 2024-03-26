import pandas as pd
import numpy as np
import plotly.figure_factory as ff
import plotly.express as px
import matplotlib.pyplot as plt
import plotly.graph_objects as go

def transform_fig_q_score_groups_by_basecaller_n_trim_levels_lengths(df):
    # Setting tresh hold for read length that is accebtable
    read_length_treshold=50

    label_less=f"length_less_{read_length_treshold}"
    label_more = f"length_more_{read_length_treshold}"

    df['count_type'] = df['len_score_sequences_read_trimmed'].apply(
        lambda x: label_less if x < read_length_treshold else label_more)

    # Removing records, of not trimmed reads
    df=df.loc[df['trim_level']!=0,:]

    # Group by the specified columns minus file_name plus the new count_type column, and count the occurrences
    grouped_counts = df.groupby(['basecaller', 'trim_level', 'trim_i_min_scores', 'count_type']).size().reset_index(
        name='count')

    # Pivot the DataFrame
    # For the new structure: Rows indexed by 'basecaller' and 'trim_i_min_scores'
    # Columns are multi-level indices composed of 'trim_level' and 'count_type'
    pivot_table = grouped_counts.pivot_table(index=['basecaller', 'trim_i_min_scores'],
                                             columns=['trim_level', 'count_type'],
                                             values='count',
                                             fill_value=0).reset_index()
    pivot_table.sort_values(by=['basecaller', 'trim_i_min_scores'], inplace=True)



    def try_convert(x):
        try:
            return int(x)
        except ValueError:
            return x  # Return the original value if it can't be converted to an integer

    # Apply the conversion function to each value in the DataFrame
    for col in pivot_table.columns:
        pivot_table[col] = pivot_table[col].apply(try_convert)

    return pivot_table

def distribution_qscore_n_nt(fq, out_format="csv"):
    """
    Transforms MyFastq.reads into long pd.DataFrame of nucleotide score per read, where columns=["read_id", 'nt', "nt_score"]
    :param fq:
    :param out_format:
    :return:
    """
    data = [[r.id, r.letter_annotations["phred_quality"]] for r in fq.reads]
    data = [[r[0], [(nt_index + 1, nt_score) for nt_index, nt_score in enumerate(r[1])]] for r in data]
    data = [[[r[0], nt[0], nt[1]] for nt in r[1]] for r in data]

    data = [x for xs in data for x in xs]

    df_data = pd.DataFrame(data, columns=["read_id", 'nt', "nt_score"])
    df_data.to_csv(f"outputs/{fq.fastq_name}.read_nt_scores.csv", index=False)
    return df_data


def trim_read_scores(read_nt_scores, score_treshold, read_lenght_treshold=20):
    """
    The goal of this function is to trim sequences based on Q score treshold
    It trimms left most and right most nucleatides that doesnt pass the score treshold
    However, it keeps nucloetides bellow the score treshold if they are between the
    left most and right most nucleotides (marked as X) that pass the score treshold

        Turn this (X is nt that pass the q sore treshold)
        [X,9,5,10,3,5,4,X]
        [10,X,5,10,3,X,4,10]
        [10,1,5,X,X,6,4,10]

        Into this
        [X,9,5,10,3,5,4,X]
        [-,X,5,10,3,X,-,-]
        [-,-,-,X,X,-,-,-]

    And discard reads that don't pass the read_lenght treshold.
    Last sequence in the example is marked as [0]

    :param read_nt_scores: list of scores from one read
    :param score_treshold: Minimum score treshold for left most and right most nucleotide
    :param read_lenght_treshold: keep reads on and above this length
    :return: list: trimmed score sequence of the read [int, int,...]
    """
    try:
        # Find indices of elements >= min_score
        read_nt_scores = np.array(read_nt_scores)
        indices = np.where(read_nt_scores >= score_treshold)[0]

        if indices.size <= 1:
            # e.g. [1,2,10,1,3,6] when score_treshold==10
            # DO NOT TRIM
            # If no nt_score is >= score_treshold (e.g. 10), do not trim - leave as it is
            # If only one nt_score is >= score_treshold (e.g. 10), do not trim - leave as it is
            # Turi būt bent dvi triminimo vietos - kiekvienam galui po vieną

            trimmed_sequence = read_nt_scores

        elif indices[0] > 0 or indices[-1] < len(read_nt_scores) - 1:
            # TRIM
            # If left most nt above score_treshhold is not first (indices[0] > 0)
            # or right most nt where score > score_treshold is not last index (indices[-1]=len(read_nt_scores)-1)

            # Use the first and last index to slice the array, keeping nt < score_treshold in the middle
            trimmed_sequence = read_nt_scores[indices[0]:indices[-1] + 1]

            """
            Since its the only case where read is trimmed, we also check if trimmed read passes the read_lenght_treshold
            If not: set trimmed sequence to [0] - nothing left from that read after trimming
            """

            if len(trimmed_sequence) <= read_lenght_treshold:
                # This sequence will return the minimum score of 0. This will mark it as "empty read" as such read will
                # have no relieble reads
                trimmed_sequence = [0]

        # DO NOT TRIM
        # If left most and right most nt that pass the treshold are first and last in the original sequence
        else:
            trimmed_sequence = read_nt_scores

    except ValueError as f:
        """
        Not sure if this one is necessary. Should find out when used in more cases
        """
        print(f"ERROR: {f}")

        print(f"score_treshold: {score_treshold}")
        print(f"read_nt_scores:{read_nt_scores}")
        print(f"trimmed_sequence {trimmed_sequence}")
        print(f"indices {indices}")

        min_nt_score_trimmed = 0

    return trimmed_sequence


def get_data_min_scores(fq, list_tresholds):
    min_trim_nt_scores = []

    for q_threshold in list_tresholds:
        # Trimming scores
        score_sequences_read_trimmed = [trim_read_scores(x, q_threshold) for x in fq.scores]

        # Getting smallest score in trimmed sequence
        trim_i_min_scores = get_min_score_group_of_trimmed_sequeces(score_sequences_read_trimmed)

        #   Counting nucleotides in score groups
        counts_trim_i_min_scores = pd.DataFrame(trim_i_min_scores, columns=["min_score"]).value_counts().reset_index()
        counts_trim_i_min_scores = df_add_fastq_metadata(counts_trim_i_min_scores, fq, q_threshold)
        min_trim_nt_scores.append(counts_trim_i_min_scores)

    df_min_trim_nt_scores = pd.concat(min_trim_nt_scores)

    return df_min_trim_nt_scores


def get_min_score_group_of_trimmed_sequeces(score_sequences_read_trimmed):
    trim_i_min_scores = [min(x) for x in score_sequences_read_trimmed]

    # Grouping qscores
    trim_i_min_scores = np.array(trim_i_min_scores)
    conditions = [(trim_i_min_scores == 0),
                  (1 <= trim_i_min_scores) & (trim_i_min_scores < 10),
                  (10 <= trim_i_min_scores) & (trim_i_min_scores < 20),
                  (20 <= trim_i_min_scores) & (trim_i_min_scores < 30),
                  (30 <= trim_i_min_scores) & (trim_i_min_scores < 40),
                  (40 <= trim_i_min_scores)]
    choices = ["Empty read", "..9", "10..19", "20..29", "30..39", "40..."]
    trim_i_min_scores = np.select(conditions, choices, "None")
    return trim_i_min_scores


def df_add_fastq_metadata(df, fq, q_threshold):
    df.loc[:, "trim_level"] = f"{q_threshold}"
    df.loc[:, "file_name"] = fq.fastq_name.replace(".fastq", "")

    df.loc[
        df["file_name"].str.contains("trimmed"),
        "adapter trimming"] = "trimmed"
    df.loc[
        ~df["file_name"].str.contains("trimmed"),
        "adapter trimming"] = "raw"

    df.loc[:, "basecaller"] = fq.basecaller
    return df


def get_data_min_scores_and_lengths(fq, list_tresholds):
    reads_trimmed_counts_lenght_min_score = []

    for q_threshold in list_tresholds:
        # Trimming scores
        score_sequences_read = [x for x in fq.scores]
        len_score_sequences_read = [len(x) for x in score_sequences_read]
        score_sequences_read_trimmed = [trim_read_scores(x, q_threshold) for x in score_sequences_read]
        len_score_sequences_read_trimmed = [len(x) for x in score_sequences_read_trimmed]

        # Getting smallest score in trimmed sequence
        trim_i_min_scores = get_min_score_group_of_trimmed_sequeces(score_sequences_read_trimmed)

        data_min_scores_lenghts = list(
            zip(len_score_sequences_read, len_score_sequences_read_trimmed, trim_i_min_scores))

        df_counts_trim_i_min_scores_lengths = pd.DataFrame.from_records(
            data=np.array(data_min_scores_lenghts),
            columns=["len_score_sequences_read", "len_score_sequences_read_trimmed", "trim_i_min_scores"])

        df_counts_trim_i_min_scores_lengths = df_counts_trim_i_min_scores_lengths.drop(
            columns=["len_score_sequences_read"])

        df_counts_trim_i_min_scores_lengths = df_counts_trim_i_min_scores_lengths.astype(
            {'len_score_sequences_read_trimmed': int})
        """
        """
        # df_counts_trim_i_min_scores_lengths=df_counts_trim_i_min_scores_lengths.value_counts(subset=["trim_i_min_scores", "len_score_sequences_read_trimmed"]).reset_index()
        # df=df_counts_trim_i_min_scores_lengths
        # fig = px.scatter(df, y="count", x="len_score_sequences_read_trimmed", color="trim_i_min_scores", marginal_y="violin",
        #                  marginal_x="box", template="simple_white")
        # fig.show()
        # exit()
        #   Counting nucleotides in score groups
        # counts_trim_i_min_scores = pd.DataFrame(trim_i_min_scores, columns=["min_score"]).value_counts().reset_index()
        counts_trim_i_min_scores = df_add_fastq_metadata(df_counts_trim_i_min_scores_lengths, fq, q_threshold)

        reads_trimmed_counts_lenght_min_score.append(counts_trim_i_min_scores)
        # read_count=counts_trim_i_min_scores.shape[0]

    df_min_trim_nt_scores = pd.concat(reads_trimmed_counts_lenght_min_score)

    # df_min_trim_nt_scores.to_csv("inputs/df_min_trim_nt_scores.tsv", sep='\t')
    # exit()

    return df_min_trim_nt_scores
