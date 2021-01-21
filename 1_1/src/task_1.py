from Bio import SeqIO


def hamming(seq1, seq2):
    assert len(seq1) == len(seq2)
    return sum(s1 != s2 for s1, s2 in zip(seq1, seq2))


seq = SeqIO.parse("../data/gattaca.fasta", "fasta")

seq1, seq2 = seq
print(hamming(seq1.seq, seq2.seq))
