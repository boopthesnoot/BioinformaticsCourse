from typing import Union
import sys
import os


sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "1_2", "src"))
from needleman import needleman_wunsch

numeric = Union[int, float]


def generate_subst_matrix(letters, match, mismatch):
    res = {}
    for i in letters:
        for j in letters:
            if i == j:
                res[(i, j)] = match
            else:
                res[(i, j)] = mismatch
    return res


def score(
    str1: str, str2: str, deletion: numeric, insertion: numeric, match: numeric, mismatch: numeric
):
    sc = [[0 for i in range(len(str2) + 1)] for j in range(2)]
    for j in range(1, len(str2) + 1):
        sc[0][j] = sc[0][j - 1] + insertion * (j - 1)
    for i in range(1, len(str1) + 1):
        sc[1][i] = sc[0][0] + deletion * (i - 1)
        for j in range(1, len(str2) + 1):
            print(i)
            if str1[i - 1] == str2[j - 1]:
                score_sub = sc[0][j - 1] + match
            else:
                score_sub = sc[0][j - 1] + mismatch * (i - 1)
            score_del = sc[0][j] + deletion * (i - 1)
            score_ins = sc[1][j - 1] + insertion * (j - 1)
            sc[1][j] = max(score_sub, score_del, score_ins)
        sc[0] = sc[1]
    return sc[1]


def hirschberg(
    str1: str, str2: str, deletion: numeric, insertion: numeric, match: numeric, mismatch: numeric
):
    res1 = []
    res2 = []
    if len(str1) == 0:
        for i in range(len(str2)):
            res1.append("-")
            res2.append(str2[i])
    elif len(str2) == 0:
        for i in range(len(str1)):
            res2.append("-")
            res1.append(str1[i])
    elif len(str1) == 1 or len(str2) == 1:
        n1, n2, junk = needleman_wunsch(
            str1, str2, generate_subst_matrix(str1 + str2, match, mismatch), deletion
        )
        res1, res2 = list(n1), list(n2)
    else:
        score1 = score(str1[: len(str1) // 2], str2, deletion, insertion, match, mismatch)
        score2 = score(
            str(reversed(str1[len(str1) // 2:])),
            str(reversed(str2)),
            deletion,
            insertion,
            match,
            mismatch,
        )
        score1.extend(score2)
        ymid = max(zip(score1, range(len(score1))))[1]
        h1 = hirschberg(str1[: len(str1) // 2], str2[:ymid], deletion, insertion, match, mismatch)
        h2 = hirschberg(str1[len(str1) // 2:], str2[ymid:], deletion, insertion, match, mismatch)
        h1[0].extend(h2[0])
        h1[1].extend(h2[1])
        res1, res2 = h1[0], h1[1]
    return res1, res2


s1, s2 = hirschberg("AGTACGCA", "TATGC", -2, -2, 2, -1)
print("".join(s1))
print("".join(s2))
