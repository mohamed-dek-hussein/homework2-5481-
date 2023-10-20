import argparse
import numpy as np

def global_alignment(seq1, seq2, match_score, mismatch_penalty, gap_penalty, unpenalized_end_gaps=False):
    """
    Perform global alignment using the Needleman-Wunsch algorithm.

    Args:
    - seq1: The first sequence to be aligned.
    - seq2: The second sequence to be aligned.
    - match_score: The score for a match.
    - mismatch_penalty: The penalty for a mismatch.
    - gap_penalty: The penalty for introducing a gap.
    - unpenalized_end_gaps: Whether to leave start and end gaps un-penalized.

    Returns:
    - Aligned sequences and the alignment score.
    """
    # Create the score matrix
    M, N = len(seq1), len(seq2)
    score_matrix = np.zeros((M + 1, N + 1), dtype=int)

    # Initialization
    for i in range(1, M + 1):
        score_matrix[i][0] = -gap_penalty * i if not unpenalized_end_gaps else 0
    for j in range(1, N + 1):
        score_matrix[0][j] = -gap_penalty * j if not unpenalized_end_gaps else 0

    # Filling the matrix
    for i in range(1, M + 1):
        for j in range(1, N + 1):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else -mismatch_penalty)
            delete = score_matrix[i-1][j] - gap_penalty
            insert = score_matrix[i][j-1] - gap_penalty
            score_matrix[i][j] = max(match, delete, insert)

    # Traceback
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = M, N
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else -mismatch_penalty):
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and score_matrix[i][j] == score_matrix[i-1][j] - gap_penalty:
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2)), score_matrix[M][N]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Needleman-Wunsch Global Sequence Alignment")
    parser.add_argument("seq1", help="First sequence")
    parser.add_argument("seq2", help="Second sequence")
    parser.add_argument("--match_score", type=int, default=1, help="Score for a match (default: 1)")
    parser.add_argument("--mismatch_penalty", type=int, default=1, help="Penalty for a mismatch (default: 1)")
    parser.add_argument("--gap_penalty", type=int, default=2, help="Penalty for a gap (default: 2)")
    parser.add_argument("--unpenalized_end_gaps", action="store_true", help="Don't penalize gaps at the start or end of alignment")
    
    args = parser.parse_args()

    aligned_seq1, aligned_seq2, score = global_alignment(args.seq1, args.seq2, args.match_score, args.mismatch_penalty, args.gap_penalty, args.unpenalized_end_gaps)
    
    print(f"Alignment score: {score}")
    print(f"Aligned Sequence 1: {aligned_seq1}")
    print(f"Aligned Sequence 2: {aligned_seq2}")
