from Bio.SubsMat import MatrixInfo as matlist

matrix = matlist.blosum62


def get_score(a, b, subst_matrix, gap):
    if a == "-" or b == "-":
        return gap
    if (a, b) not in subst_matrix:
        return subst_matrix[(b, a)]
    else:
        return subst_matrix[(a, b)]


def needleman_wunsch(seq1, seq2, subst_matrix, gap):
    n = len(seq1)
    m = len(seq2)
    score = [[0 for _ in range(n + 1)] for j in range(m + 1)]
    for i in range(m + 1):
        score[i][0] = gap * i
    for i in range(n + 1):
        score[0][i] = gap * i
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + get_score(seq1[j - 1], seq2[i - 1], subst_matrix, gap)
            insertion = score[i - 1][j] + gap
            deletion = score[i][j - 1] + gap
            score[i][j] = max(match, insertion, deletion)
    seq1_a = []
    seq2_a = []
    i = m
    j = n

    while i > 0 and j > 0:
        cur_score = score[i][j]
        diag_score = score[i - 1][j - 1]
        up_score = score[i][j - 1]

        if cur_score == diag_score + get_score(seq1[j - 1], seq2[i - 1], subst_matrix, gap):
            seq1_a.append(seq1[j - 1])
            seq2_a.append(seq2[i - 1])
            i -= 1
            j -= 1
        elif cur_score == up_score + gap:
            seq1_a.append(seq1[j - 1])
            seq2_a.append("-")
            j -= 1
        else:
            seq1_a.append("-")
            seq2_a.append(seq2[i - 1])
            i -= 1
    while i > 0:
        seq2_a.append(seq2[i - 1])
        seq1_a.append("-")
        i -= 1
    while j > 0:
        seq2_a.append("-")
        seq1_a.append(seq1[j - 1])
        j -= 1
    return "".join(reversed(seq1_a)), "".join(reversed(seq2_a)), score[-1][-1]


for i in needleman_wunsch("ATTACA", "ATGCT", matrix, -3):
    print(i)

