import numpy as np


def local_alignment(s, t, match, mismatch, gap):
    """
    Local alignment of two strings s and t.
    :param s: string 1
    :param t: string 2
    :param match: match score
    :param mismatch: mismatch score
    :param gap: gap score
    :return: dp_matrix, score, and list of aligned strings
    """
    # length of s and t
    m = len(s)
    n = len(t)

    # initialize dp_matrix
    dp = [[0 for i in range(m + 1)] for j in range(n + 1)]

    for i in range(m + 1):
        dp[0][i] = 0

    for i in range(n + 1):
        dp[i][0] = 0

    # fill in dp_matrix
    mv = 0
    mx, my = 0, 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if s[j - 1] == t[i - 1]:
                dp[i][j] = max(0, dp[i - 1][j - 1] + match, dp[i - 1][j] + gap, dp[i][j - 1] + gap)
            else:
                dp[i][j] = max(0, dp[i - 1][j - 1] + mismatch, dp[i - 1][j] + gap, dp[i][j - 1] + gap)
            if dp[i][j] > mv:
                mx, my = i, j
                mv = dp[i][j]


    print(np.array(dp))
    print('score : ', mv)

    l = []
    # find the path
    def r(i, j, a1, a2):
        # print(a1, a2)
        if dp[i][j] == 0:
            l.append((a1, a2))
            return
        # if the score is equal to the score in the diagonal cell+match, then move diagonally and the characters are
        # the same
        if i > 0 and j > 0 and s[j - 1] == t[i - 1] and dp[i][j] == dp[i - 1][j - 1] + match:
            r(i - 1, j - 1, s[j - 1] + a1, t[i - 1] + a2)
        # if the score is equal to the score in the diagonal cell+mismatch, then move up and the characters are
        # different
        if i > 0 and j > 0 and s[j - 1] != t[i - 1] and dp[i][j] == dp[i - 1][j - 1] + mismatch:
            r(i - 1, j - 1, s[j - 1] + a1, t[i - 1] + a2)
        # if the score is equal to the score in the up cell+gap, then move up and add gap to a1
        if i > 0 and dp[i][j] == dp[i - 1][j] + gap:
            r(i - 1, j, '-' + a1, t[i - 1] + a2)
        # if the score is equal to the score in the left cell+gap, then move left and add gap to a2
        if j > 0 and dp[i][j] == dp[i][j - 1] + gap:
            r(i, j - 1, s[j - 1] + a1, '-' + a2)

    r(mx, my, '', '')
    print('optimal alignments:', l, '\n')
    return dp, dp[mx][my], l


if __name__ == '__main__':
    local_alignment('ATCAGAGTA', 'TTCAGTA', 2, -1, -1)
