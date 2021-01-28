from Bio import SeqIO


def levenshtein(seq1, seq2):
    size1 = len(seq1)
    size2 = len(seq2)
    m = [[0] * size2 for _ in range(size1)]
    for i in range(size1):
        m[i][0] = i
    for i in range(size2):
        m[0][i] = i

    for i in range(1, size1):
        for j in range(1, size2):
            if seq1[i - 1] == seq2[j - 1]:
                m[i][j] = min(m[i - 1][j] + 1, m[i - 1][j - 1], m[i][j - 1] + 1)
            else:
                m[i][j] = min(m[i - 1][j] + 1, m[i - 1][j - 1] + 1, m[i][j - 1] + 1)
    return m[-1][-1]


seq = SeqIO.parse("../data/gattaca.fasta", "fasta")

seq1, seq2 = seq
print(levenshtein(seq1.seq, seq2.seq))
