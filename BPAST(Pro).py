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

# 定义位置判定函数 (Define the position determination function)
def check_position_relation(start1, end1, start2, end2):
    if end1 < start2:
        return "B"  # 子序列1在子序列2之前 (Subsequence 1 before subsequence 2)
    elif start1 > end2:
        return "A"  # 子序列1在子序列2之后 (Subsequence 1 after subsequence 2)
    else:
        return "E"  # 子序列1与子序列2相交或包含关系 (Subsequence 1 intersects or contains a relationship with subsequence 2)

# 定义比较序列的函数 (A function that defines a comparison sequence)
def compare_sequences(m_seqs1, m_seqs2, n_file, positions1, positions2, acids1, acids2, k1, k2, S1, S2, percent_range1, percent_range2, position_comparison, output_file):
    results = []
    with open(n_file, "r") as f:
        n_sequences = [seq_data.split('\n', 1) for seq_data in f.read().split('>')[1:]]
    n_sequences = [(header.split()[0], seq.replace('\n', '')) for header, seq in n_sequences]

    def process_comparison(m_seqs, percent_range, k, positions, acids, S):
        match_results = []
        for m_seq in m_seqs:
            for header, n_seq in n_sequences:
                n_seq_length = len(n_seq)
                subseq_start, subseq_end = L(n_seq_length, percent_range)
                lengths = [len(m_seq) + i for i in range(-k, k + 1)]
                best_match = None

                for length in lengths:
                    subsequences = generate_subsequences(n_seq[subseq_start:subseq_end], length)
                    for start, subseq in enumerate(subsequences, subseq_start):
                        score = compute_score(m_seq, subseq)
                        if best_match is None or score > best_match['score']:
                            best_match = {
                                'header': header,
                                'n_seq_length': n_seq_length,
                                'start': start + 1,
                                'end': start + length,
                                'subseq': subseq,
                                'score': score,
                                'percentage_score': (score / len(m_seq)) * 100
                            }

                match_count = sum(
                    (pos > 0 and best_match['subseq'][pos - 1] == ac) or
                    (pos <= 0 and best_match['subseq'][pos] == ac)
                    for pos, ac in zip(positions, acids)
                ) if best_match else 0

                if match_count >= S:
                    match_results.append([
                        m_seq,
                        best_match['header'],
                        best_match['n_seq_length'],
                        best_match['start'],
                        best_match['end'],
                        best_match['subseq'],
                        best_match['percentage_score']
                    ])
        return match_results

    # 处理第一组参数 (Processing a first set of parameters)
    results1 = process_comparison(m_seqs1, percent_range1, k1, positions1, acids1, S1)

    # 处理第二组参数 (Processing the second set of parameters)
    results2 = process_comparison(m_seqs2, percent_range2, k2, positions2, acids2, S2)

    # 合并两个结果并进行位置判定 (Combine the two results and make position determination)
    for res1 in results1:
        for res2 in results2:
            if res1[1] == res2[1]:  # 确保来自同一条序列 (Make sure they come from the same sequence)
                position_relation = check_position_relation(res1[3], res1[4], res2[3], res2[4])
                if position_relation in position_comparison:
                    combined_result = [
                        res1[1],  # Sequence ID
                        res1[2],  # Length
                        res1[0],  # Query m_seq1
                        res1[3],  # Start1
                        res1[4],  # End1
                        res1[5],  # Matching Sub-sequence 1
                        res1[6],  # Score1 (%)
                        res2[0],  # Query m_seq2
                        res2[3],  # Start2
                        res2[4],  # End2
                        res2[5],  # Matching Sub-sequence 2
                        res2[6],  # Score2 (%)
                        position_relation  # Relative Position
                    ]
                    results.append(combined_result)

    results_df = pd.DataFrame(results, columns=[
        "Sequence ID", "Length", "Query m_seq1", "Start1", "End1", "Matching Sub-sequence 1", "Score1 (%)",
        "Query m_seq2", "Start2", "End2", "Matching Sub-sequence 2", "Score2 (%)", "Position Relation"
    ])

    # 创建保存参数的DataFrame (Create a DataFrame to save parameters)
    parameters_data = {
        'Parameter': ['gap_penalty', 'match_score', 'mismatch_penalty', 'k1', 'k2', 'S1', 'S2',
                      'percent_range1', 'percent_range2', 'positions1', 'positions2', 'acids1', 'acids2',
                      'm_seqs1', 'm_seqs2', 'n_file', 'output_file', 'position_comparison'],
        'Value': [
            gap_penalty,
            match_score,
            mismatch_penalty,
            k1,
            k2,
            S1,
            S2,
            str(percent_range1),
            str(percent_range2),
            positions1,
            positions2,
            acids1,
            acids2,
            ' '.join(m_seqs1),
            ' '.join(m_seqs2),
            n_file,
            output_file,
            position_comparison
        ]
    }
    parameters_df = pd.DataFrame(parameters_data)

    # 写入结果到Excel文件的不同sheets (Write the results to different sheets of Excel file)
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        results_df.to_excel(writer, index=False, sheet_name='Combined Matches')
        parameters_df.to_excel(writer, index=False, sheet_name='Parameters')

    print(f'Results and parameters have been saved to: {output_file}')
    return output_file



# 参数设置 (Parameter settings)
k1 = 0
S1 = 3
percent_range1 = [43, 72]
positions1 = [1,2,3]
acids1 = ['V', 'W','D']
m_seqs1 = ['VWDQ', 'VWDH', 'VWDK']

k2 = 0
S2 = 2
percent_range2 = [58, 100]
positions2 = [5, 7]
acids2 = ['S', 'S']
m_seqs2 = ['EVGGSCSPHAHGR', 'EVDGSCSRRAPGR', 'GASGSNSGRAPSC', 'RVLPSASRRGPGQ', 'GVGASNSGHSPGA',
           'AVGGSDSVRAHSK', 'AVDHSYSPRAGQR', 'DVGGSNSRRAGQG', 'EVRGSTSGSAGQG']

position_comparison = ["B"]  # 可选值: "B", "A", "E" (Optional values: "B", "A", "E a" and "e")

n_file = r"C:\Users\27404\Arabidopsis_protein.fa"
output_file = r"C:\Users\27404\Arabidopsis_protein_SCOOP.xls"

# 执行函数 (Executive function)
result_file_path = compare_sequences(m_seqs1, m_seqs2, n_file, positions1, positions2, acids1, acids2, k1, k2, S1, S2, percent_range1, percent_range2, position_comparison, output_file)