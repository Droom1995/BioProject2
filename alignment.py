import numpy as np


def needlman_wunsch(seq1, seq2, letters, blosum62, gap_p):
    n = len(seq1)
    m = len(seq2)
    dp = np.zeros((m + 1, n + 1))
    path = np.zeros((m + 1, n + 1))
    dp[0][0] = 0
    for i in range(1, n + 1):
        dp[0][i] = gap_p * i
    for j in range(1, m + 1):
        dp[j][0] = gap_p * j
    for j in range(1, n + 1):
        for i in range(1, m + 1):
            match = dp[i - 1][j - 1] + match_score(seq2[i - 1], seq1[j - 1], letters, blosum62)
            delete = dp[i - 1][j] + gap_p
            insert = dp[i][j - 1] + gap_p
            dp[i][j] = max(match, delete, insert)
            if dp[i][j] == delete:
                path[i][j] = 1  # up
            if dp[i][j] == insert:
                path[i][j] = 2  # left
            if dp[i][j] == match:
                path[i][j] = 3  # diagonal
    align1, align2 = find_alignment(seq1=seq2, seq2=seq1, path=path, max_i=m, max_j=n, global_sym='-')

    return align1, align2, int(dp[m][n])


def smith_waterman(seq1, seq2, letters, blosum62, gap_penalty):
    m, n = len(seq1), len(seq2)

    dp = np.zeros((m + 1, n + 1))
    path = np.zeros((m + 1, n + 1))

    max_score = 0
    max_i, max_j = 0, 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_diagonal = dp[i - 1][j - 1] + match_score(seq1[i - 1], seq2[j - 1], letters, blosum62)
            score_up = dp[i][j - 1] + gap_penalty
            score_left = dp[i - 1][j] + gap_penalty
            dp[i][j] = max(0, score_left, score_up, score_diagonal)
            if dp[i][j] == 0:
                path[i][j] = 0 # cut path
            if dp[i][j] == score_left:
                path[i][j] = 1  # up
            if dp[i][j] == score_up:
                path[i][j] = 2  # left
            if dp[i][j] == score_diagonal:
                path[i][j] = 3  # diagonal
            if dp[i][j] >= max_score:
                max_i = i
                max_j = j
                max_score = dp[i][j]
    align1, align2 = find_alignment(seq1, seq2, path, max_i, max_j)
    return align1, align2, int(dp[max_i][max_j])


def match_score(a, b, letters, blosum62):
    return blosum62[letters.index(a)][letters.index(b)]


def find_alignment(seq1, seq2, path, max_i, max_j,global_sym=''):
    align1, align2 = '', ''
    i, j = max_i, max_j
    while path[i][j] != 0:
        if path[i][j] == 3:
            align1 += seq1[i - 1] if i > 0 else global_sym
            align2 += seq2[j - 1] if j > 0 else global_sym
            i -= 1
            j -= 1
        elif path[i][j] == 2:
            align1 += '-'
            align2 += seq2[j - 1]
            j -= 1
        elif path[i][j] == 1:
            align1 += seq1[i - 1]
            align2 += '-'
            i -= 1
    align1 = align1[::-1]
    align2 = align2[::-1]
    return align1, align2
