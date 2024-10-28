import itertools
from difflib import SequenceMatcher
import pandas as pd

# 设置得分参数 (Set scoring parameters)
gap_penalty = -1
match_score = 1
mismatch_penalty = 0


# 定义函数来计算氨基酸得分 (Define a function to calculate the amino acid score)
def get_amino_acid_score(a, b):
    if a == b:
        return match_score
    elif a == '-' or b == '-':
        return gap_penalty
    else:
        return mismatch_penalty


# 定义函数来计算两个序列之间的总得分 (Define a function to calculate the total score between two sequences)
def compute_score(seq1, seq2):
    score = 0
    for a, b in zip(seq1, seq2):
        score += get_amino_acid_score(a, b)
    return score


# 定义生成可能子序列的函数 (Defines the function that generates the possible subsequences)
def generate_subsequences(sequence, length):
    subsequences = []
    for i in range(len(sequence) - length + 1):
        subsequences.append(sequence[i:i + length])
    return subsequences


# 定义函数来计算相似度 (Define a function to calculate the similarity)
def similarity(seq1, seq2):
    return SequenceMatcher(None, seq1, seq2).ratio()


def L(sequence_length, percent_range):
    start = int(sequence_length * percent_range[0] / 100)
    end = int(sequence_length * percent_range[1] / 100)
    return start, end


# 定义检查子序列是否符合条件的函数 (Defines a function that checks whether a subsequence meets a condition)
def check_subseq_conditions(subseq, positions, acids, S):
    match_count = sum(
        (pos > 0 and subseq[pos - 1] == ac) or
        (pos <= 0 and subseq[pos] == ac)
        for pos, ac in zip(positions, acids)
    )
    return match_count >= S


# 定义比较序列的函数 (A function that defines a comparison sequence)
def compare_sequences(m_seqs, n_file, positions1, acids1, k1, S1, positions2, acids2, k2, S2, positions3, acids3, k3,
                      S3, percent_range, output_file):
    results = []
    with open(n_file, "r") as f:
        n_sequences = [seq_data.split('\n', 1) for seq_data in f.read().split('>')[1:]]
    n_sequences = [(header.split()[0], seq.replace('\n', '')) for header, seq in n_sequences]

    for m_seq in m_seqs:
        for header, n_seq in n_sequences:
            n_seq_length = len(n_seq)
            subseq_start, subseq_end = L(n_seq_length, percent_range)
            lengths = [len(m_seq) + i for i in range(-max(k1, k2, k3), max(k1, k2, k3) + 1)]

            for length in lengths:
                subsequences = generate_subsequences(n_seq[subseq_start:subseq_end], length)
                for start, subseq in enumerate(subsequences, subseq_start):
                    score = compute_score(m_seq, subseq)
                    percentage_score = (score / len(m_seq)) * 100

                    # 检查三个条件 (Check three conditions)
                    condition1_met = check_subseq_conditions(subseq, positions1, acids1, S1)
                    condition2_met = check_subseq_conditions(subseq, positions2, acids2, S2)
                    condition3_met = check_subseq_conditions(subseq, positions3, acids3, S3)

                    if condition1_met and condition2_met and condition3_met:
                        results.append([
                            m_seq,
                            header,
                            n_seq_length,
                            start + 1,
                            start + length,
                            subseq,
                            percentage_score
                        ])

    results_df = pd.DataFrame(results, columns=[
        "Query m_seq", "Sequence ID", "Length", "Start", "End",
        "Matching Sub-sequence", "Score (%)"
    ])
    results_df_sorted = results_df.sort_values(by="Score (%)", ascending=False)
    results_df_unique = results_df_sorted.drop_duplicates(subset="Sequence ID", keep="first")

    # 创建保存参数的DataFrame (Create a DataFrame to save parameters)
    parameters_data = {
        'Parameter': ['gap_penalty', 'match_score', 'mismatch_penalty', 'k1', 'S1', 'k2', 'S2', 'k3', 'S3',
                      'percent_range', 'positions1', 'acids1', 'positions2', 'acids2', 'positions3', 'acids3', 'm_seqs',
                      'n_file', 'output_file'],
        'Value': [
            gap_penalty, match_score, mismatch_penalty, k1, S1, k2, S2, k3, S3, str(percent_range),
            positions1, acids1, positions2, acids2, positions3, acids3, ' '.join(m_seqs), n_file, output_file
        ]
    }
    parameters_df = pd.DataFrame(parameters_data)

    # 计算序列之间的相似度 (Calculate the similarity between sequences)
    similarities = [(seq1, seq2, similarity(seq1, seq2)) for seq1, seq2 in itertools.combinations(m_seqs, 2)]
    similarities_df = pd.DataFrame(similarities, columns=['Sequence 1', 'Sequence 2', 'Similarity'])

    # 写入结果到Excel文件的不同sheets (Write the results to different sheets of Excel file)
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        results_df_sorted.to_excel(writer, index=False, sheet_name='All Matches')
        results_df_unique.to_excel(writer, index=False, sheet_name='Unique Best Matches')
        parameters_df.to_excel(writer, index=False, sheet_name='Parameters')
        similarities_df.to_excel(writer, index=False, sheet_name='Similarity Analysis')

    print(f'Results and parameters have been saved to: {output_file}')
    return output_file


# 参数设置 (Parameter settings)
k1 = 0
S1 = 15
positions1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
acids1 = ['A','C','C','A','C','A','G','C','A','G','C','C','T','C','T','G','A','A','G','T']
k2 = 0
S2 = 2
positions2 = [22, 23]
acids2 = ['G', 'G']
k3 = 0
S3 = 0
positions3 = [21]
acids3 = ['A']
percent_range = [0, 100]
m_seqs = ['ACCACAGCAGCCTCTGAAGTAGG']
n_file = r"C:\Users\27404\Arabidopsis_Genome.fa"
output_file = r"C:\Users\27404\Arabidopsis_CRIPSR_Target.xlsx"

# 执行函数 (Executive function)
result_file_path = compare_sequences(m_seqs, n_file, positions1, acids1, k1, S1, positions2, acids2, k2, S2, positions3,
                                     acids3, k3, S3, percent_range, output_file)