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

# 定义比较序列的函数 (A function that defines a comparison sequence)
def compare_sequences(m_seqs, n_file, positions, acids, k, S, output_file):
    results = []
    with open(n_file, "r") as f:
        n_sequences = [seq_data.split('\n', 1) for seq_data in f.read().split('>')[1:]]
    n_sequences = [(header.split()[0], seq.replace('\n', '')) for header, seq in n_sequences]

    for m_seq in m_seqs:
        for header, n_seq in n_sequences:
            n_seq_length = len(n_seq)  # 获取 n 序列的长度
            subseq_start, subseq_end = L(n_seq_length, percent_range)  # 获得比对位置的起始和结束 (Get the start and end of the comparison position)
            lengths = [len(m_seq) + i for i in range(-k, k + 1)]
            best_match = None

            for length in lengths:
                subsequences = generate_subsequences(n_seq[subseq_start:subseq_end], length)
                for start, subseq in enumerate(subsequences, subseq_start):# 从 subseq_start 开始枚举 (Enumeration begins with subseq_start)
                    score = compute_score(m_seq, subseq)
                    if best_match is None or score > best_match['score']:
                        best_match = {
                            'header': header,
                            'n_seq_length': n_seq_length,  # 保存 n 序列的长度 (Save the length of the n sequence)
                            'start': start + 1,
                            'end': start + length,
                            'subseq': subseq,
                            'score': score,
                            'percentage_score': (score / len(m_seq)) * 100
                        }

            # 修改匹配共识位置的逻辑以考虑新的S值 (Modify the logic that matches the consensus position to take into account the new S value)
            match_count = sum(
                (pos > 0 and best_match['subseq'][pos - 1] == ac) or
                (pos <= 0 and best_match['subseq'][pos] == ac)
                for pos, ac in zip(positions, acids)
            ) if best_match else 0

            if match_count >= S:  # 添加 S 值的匹配条件 (Add matching condition of s value)
                results.append([
                    m_seq,
                    best_match['header'],
                    best_match['n_seq_length'],  # 添加 n 序列长度到结果 (Add n sequence length to the result)
                    best_match['start'],
                    best_match['end'],
                    best_match['subseq'],
                    best_match['percentage_score']
                ])

            match_consensus_positions = all(
                (pos > 0 and best_match['subseq'][pos - 1] == ac) or
                (pos <= 0 and best_match['subseq'][pos] == ac)
                for pos, ac in zip(positions, acids)
            ) if best_match else False

            if match_consensus_positions:
                results.append([
                    m_seq,
                    best_match['header'],
                    best_match['start'],
                    best_match['end'],
                    best_match['subseq'],
                    best_match['percentage_score']
                ])

    results_df = pd.DataFrame(results, columns=[
        "Query m_seq", "Sequence ID", "Length", "Start", "End",
        "Matching Sub-sequence", "Score (%)"
    ])
    results_df_sorted = results_df.sort_values(by="Score (%)", ascending=False)
    results_df_unique = results_df_sorted.drop_duplicates(subset="Sequence ID", keep="first")

    # 创建保存参数的DataFrame (Create a DataFrame to save parameters)
    parameters_data = {
        'Parameter': ['gap_penalty', 'match_score', 'mismatch_penalty', 'k', 'S','percent_range','positions', 'acids', 'm_seqs', 'n_file',
                      'output_file'],
        'Value': [
            gap_penalty,
            match_score,
            mismatch_penalty,
            k,
            S,
            str(percent_range),
            positions,
            acids,
            ' '.join(m_seqs),  # 将序列列表转换为字符串 (Convert a sequence list to a string)
            n_file,
            output_file
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
k = 0
S = 4
percent_range = [0, 100]
positions = [ 4, 5, -2, -1]
acids = ['S','G','G','H']
m_seqs = ['RLASGPSPRGRGH','VKHSGPSPSGPGH','GKHSGPSTSGPGH','KLASGPSRRGCGH','TMASGPSRRGAGH','ILASGPNKRGRGH','RLASGPSRKGRGH','RLASGPSRRGRGH','KLVSGPSRSSCGH','RLASGSSRRGRGH']
n_file = r"C:\Yxs-Doctor work\TAIR10.0.fa"
output_file = r"C:\Yxs-Doctor work\PIP.xls"


# 执行函数 (Executive function)
result_file_path = compare_sequences(m_seqs, n_file, positions, acids, k, S, output_file)
