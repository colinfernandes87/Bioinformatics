""" THE UNIVERSAL TOOL FOR ALIGNMENT """
""" Author: COLIN MANUEL FERNANDES """
""" SUBJECT: BIOINFORMATICS """
""" 914777 """

# Libraries
import re
import numpy as np
import sys

""" GLOBAL ALIGNMENT OF THE GIVEN SEQUENCES """
""" NEEDLEMAN - WUNSCH ALGORITHM """
""" APPLICATIONS: Compare two genes of same sequences """


# --------------------------------------------------------------------

def globl(seq1, seq2):
    m, n = len(seq1), len(seq2)  # length of two sequences

    # Generate DP table and traceback path pointer matrix
    score = zeros((m + 1, n + 1))  # the DP table

    # Calculate DP table
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i
    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + match_score(seq1[i - 1], seq2[j - 1])
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert)

    # Traceback and compute the alignment
    align1, align2 = '', ''
    i, j = m, n  # start from the bottom right cell
    while i > 0 and j > 0:  # end crossing the top or the left edge
        score_current = score[i][j]
        score_diagonal = score[i - 1][j - 1]
        score_up = score[i][j - 1]
        score_left = score[i - 1][j]

        if score_current == score_diagonal + match_score(seq1[i - 1], seq2[j - 1]):
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += seq1[i - 1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + gap_penalty:
            align1 += '-'
            align2 += seq2[j - 1]
            j -= 1

    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq1[i - 1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j - 1]
        j -= 1

    print("\n\nBelow is the Global Alignment\n")
    finalize(align1, align2, score)


# --------------------------------------------------------------------


""" LOCAL ALIGNMENT OF THE GIVEN SEQUENCES """
""" SMITH WALTERMAN ALGORITHM """
""" APPLICATIONS:Searching local similarities """


def locl(seq1, seq2):
    m, n = len(seq1), len(seq2)  # length of two sequences
    max_score = 0  # initial maximum score in DP table

    # Generate DP table and traceback path pointer matrix
    score = zeros((m + 1, n + 1))  # the DP table

    # Calculate DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + match_score(seq1[i - 1], seq2[j - 1])
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(0, match, delete, insert)
            if score[i][j] > max_score:
                max_i = i
                max_j = j
                max_score = score[i][j]

    align1, align2 = '', ''  # initial sequences

    i, j = max_i, max_j  # indices of path starting point

    # Traceback and compute the alignment
    while score[i][j] != 0:  # end crossing the top or the left edge
        score_current = score[i][j]
        score_diagonal = score[i - 1][j - 1]
        score_up = score[i][j - 1]
        score_left = score[i - 1][j]

        if score_current == score_diagonal + match_score(seq1[i - 1], seq2[j - 1]):
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += seq1[i - 1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + gap_penalty:
            align1 += '-'
            align2 += seq2[j - 1]
            j -= 1

    print("\n\nBelow is the Local Alignment\n")
    finalize(align1, align2, score)


# --------------------------------------------------------------------


""" SEMI-GLOBAL ALIGNMENT OF THE GIVEN SEQUENCES """
""" FREE END GAP MODIFICATION OF NEEDLEMAN - WUNSCH ALGORITHM """
""" APPLICATIONS: sequence assembly, overlap detection """


def semiglobal_alignment(seq1, seq2):
    m = len(seq1) + 1
    n = len(seq2) + 1

    # Initialization
    # scores matrix
    score = [[0 for i in range(m)] for j in range(n)]
    # traceback matrix
    traceback = [[0 for i in range(m)] for j in range(n)]  # to store the trace back path

    # Dynamic programing table
    for i in range(1, n):
        for j in range(1, m):
            gap_penalty_i, gap_penalty_j = gap_penalty, gap_penalty
            if i == m - 1:
                gap_penalty_i = 0
            if j == n - 1:
                gap_penalty_j = 0

            if seq1[j - 1] == seq2[i - 1]:
                match = score[i - 1][j - 1] + match_award
            else:
                match = score[i - 1][j - 1] + mismatch_penalty
            score_left = score[i][j - 1] + gap_penalty_i
            score_above = score[i - 1][j] + gap_penalty_j
            score[i][j] = max(score_left, score_above, match)

            if score[i][j] == match:
                traceback[i][j] = 1  # 1 means trace diagonally
            elif score[i][j] == score_left:
                traceback[i][j] = 2  # 2 means trace to the left
            else:
                traceback[i][j] = 3  # 3 means trace to the top
    # Initialization
    align1, align2 = '', ''

    # stores the number of column the max number
    max_j = score[-1].index(max(score[-1]))

    """ ----COMPUTE THE TRACEBACK ALIGNMENT---- """
    while i >= 0 and j >= 0:
        if j > max_j:
            align1 = seq1[j - 1] + align1
            align2 = '-' + align2
            j -= 1
            continue

        if traceback[i][j] == 1:  # 1 means to trace diagonally
            align1 = seq1[j - 1] + align1
            align2 = seq2[i - 1] + align2
            i -= 1
            j -= 1
        elif traceback[i][j] == 2:  # 2 means trace to the left
            align1 = seq1[j - 1] + align1
            align2 = '-' + align2
            j -= 1
        elif traceback[i][j] == 0:  # 0 means trace to the left
            if j == 0:
                align1 = '-' + align1
                align2 = seq2[i - 1] + align2
                j -= 1
            else:
                align1 = seq1[j - 1] + align1
                align2 = '-' + align2
                i -= 1
        else:  # 3 means trace to the top
            align1 = '-' + align1
            align2 = seq2[i - 1] + align2
            i -= 1

    print("\n\nBelow is the SemiGlobal Alignment\n")
    print("s: ", align1)
    print("t: ", align2)
    print("Dynamic Programing Matrix: \n", np.matrix(score))
    print("\nMaximum Score: ", np.matrix(score).max())


def zeros (shape):
    retval = []
    for x in range(shape[0]):
        retval.append([])
        for y in range(shape[1]):
            retval[-1].append(0)
    return retval


def match_score (alpha, beta):
    if alpha == beta:
        return match_award
    elif alpha == '-' or beta == '-':
        return gap_penalty
    else:
        return mismatch_penalty


def finalize (align1, align2, scores):
    align1 = align1[::-1]  # reverse sequence 1
    align2 = align2[::-1]  # reverse sequence 2

    # calcuate identity, score and aligned sequeces
    score = 0
    for i in range(0, len(align1)):
        # if two AAs are the same, then output the letter
        if align1[i] == align2[i]:
            score += match_score(align1[i], align2[i])

        # if they are not identical and none of them is gap
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-':
            score += match_score(align1[i], align2[i])

        # if one of them is a gap, output a space
        elif align1[i] == '-' or align2[i] == '-':
            score += gap_penalty

    print("s: ", align1)
    print("t: ", align2)
    print("\nDynamic Programing Matrix: \n", np.matrix(scores))
    print("\nAlignment Score: ", score)


" USER OPTIONS"


def retest():
    while True:
        action = input("Do you choose another test 'Y' or 'N' ").upper()
        if action == "N":
            print("Thank you for Using UniAlign \n "
                  "Please contact the developer (colinmanuel.fernandes@studenti.unimi.it) for issues\n")
            sys.exit()
        else:
            print("Welcome back!!\n")
            main()


" VALIDATION OF USER INPUT"


def validate_char(ipt):
    if not re.match("^[A-Z]*$", ipt):
        print("Error! Only letters a-z allowed!")
        sys.exit()


" VALIDATION OF USER INPUT"


def validate_int(ipt):
    if not re.match("^[-+]?([1-9]\d*|0)$", ipt):
        print("Error! Only integers with '+' or '-' allowed!")
        sys.exit()


# --------------------------------------------------------------------


" LONGEST COMMON SUB- SEQUENCE"


def lcs(seq1, seq2):
    score = [[0 for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]
    # Initialization
    for i, x in enumerate(seq1):
        for j, y in enumerate(seq2):
            if x == y:
                score[i + 1][j + 1] = score[i][j] + 1
            else:
                score[i + 1][j + 1] = max(score[i + 1][j], score[i][j + 1])
    # read the LCS out from the matrix
    seq = ""
    x, y = len(seq1), len(seq2)
    while x != 0 and y != 0:
        if score[x][y] == score[x - 1][y]:
            x -= 1
        elif score[x][y] == score[x][y - 1]:
            y -= 1
        else:
            assert seq1[x - 1] == seq2[y - 1]
            seq = seq1[x - 1] + seq
            x -= 1
            y -= 1
    print("Longest Sub-Sequence of the string : ", seq)
    print("Length of the Sub-Sequence :", len(seq))
    print("Edit distance required to change 1st string to 2nd  :", len(seq1)+ len(seq2) - 2*len(seq))

# --------------------------------------------------------------------

" PROGRAM MAIN FOR USER INPUT"


def main():
    while True:
        action = input(
            "Choose your alignment :\n"
            "'G' for Needleman-Wunsch (Global) \n"
            "'L' for Smith-Waterman (Local)\n"
            "'S' for End gap Modified(Semi Global)\n"
            "'C' for Longest Common Sub Sequence\n"
            "'A' for all Alignments\n"
            "'E' Exit? \n").upper()
        if action not in "GLSACE" or len(action) != 1:
            print("This is new for me and I don't know how to do that yet,"
                  "\n you can contact the developer (colinmanuel.fernandes@studenti.unimi.it) for implementation to "
                  "add a new algorithm")
        if action == 'E':
            sys.exit()
        sequence1 = input("Enter your 1st Sequence s:").upper()
        validate_char(sequence1)
        sequence2 = input("Enter your 2nd Sequence t:").upper()
        validate_char(sequence2)
        if action != 'C':
            global match_award
            global mismatch_penalty
            global gap_penalty
            match_award = int(input("Input Match Award"))
            # validate_int(match_award)
            mismatch_penalty = int(input("Input miss match Penalty "))
            # validate_int(mismatch_penalty)
            gap_penalty = int(input("Input gap penalty"))
            # validate_int(gap_penalty)
        if action == 'S':
            semiglobal_alignment(sequence1, sequence2)
        elif action == 'G':
            globl(sequence1, sequence2)
        elif action == 'L':
            locl(sequence1, sequence2)
        elif action == 'C':
            lcs(sequence1, sequence2)
        elif action == 'A':
            globl(sequence1, sequence2)
            locl(sequence1, sequence2)
            semiglobal_alignment(sequence1, sequence2)
            lcs(sequence1, sequence2)
        retest()


" EXECUTION POINT "
if __name__ == "__main__":
    print("WELCOME TO UNI-ALIGN: UNIVERSAL TOOL FOR GLOBAL ,LOCAL , SEMIGLOBAL ALIGNMENT ")
    main()
    sys.exit()
