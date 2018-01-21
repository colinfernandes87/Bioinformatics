# Bioinformatics

In bioinformatics, a sequence alignment is a way of arranging the sequences of DNA, RNA, or protein to identify regions of similarity that may be a consequence of functional, structural, or evolutionary relationships between the sequences.[1] Aligned sequences of nucleotide or amino acid residues are typically represented as rows within a matrix. Gaps are inserted between the residues so that identical or similar characters are aligned in successive columns. Sequence alignments are also used for non-biological sequences, such as calculating the edit distance cost between strings in a natural language or in financial data


Implementation of algorithms for Bio Informatics 

Uni-Align

Description: A General application for checking all the sequence alignments using dynamic programing and also to compute the edit distance.

The application is developed in Python and consists of following implementations
1)Global Alignment
2)Local Alignment
3)Semi Global (Glocal) Alignment
4)Longest Common Subsequence (find the edit distance between the strings)

Below is the Theoretical description of each of the components in the code

1. Global Alignment: 
Definition1: A (global) alignment of two strings S and T is obtained by first inserting
chosen spaces (or dashes), either into or at the ends of S and T, and then placing the
two resulting strings one above the other so that every character or space in either string
is opposite a unique character or a unique space in the other string.

Algorithm Used: Needleman–Wunsch algorithm2 based on Dynamic Programing

2. Local Alignment: 
Definition3: A local alignment between S and T is an alignment between a substring of S and a substring
of T. In this section we present an algorithm to find the highest scoring local alignments
between two sequences.

Algorithm Used: Smith–Waterman algorithm based on Dynamic Programing

3. Semiglobal Alignment: 
Definition4: In a semiglobal comparison, we score alignments ignoring some of the end spaces in the
sequences. An interesting characteristic of the basic dynamic programming algorithm is
that we can control the penalty associated with end spaces by doing very simple modifications
to the original scheme.

Algorithm Used: Free end gap modification of Needleman–Wunsch Algorithm based on Dynamic Programing

4. Longest Common Subsequence: 
Definition5: Given two sequences X and Y, we say that a sequence Z is a common subsequence
of X and Y if Z is a subsequence of both X and Y . For example, if
X = _A, B, C, B, D, A, B_ and Y = _B, D, C, A, B, A_, the sequence B, C, A
is a common subsequence of both X and Y. In the longest-common-subsequence problem, we are given two sequences X =_x1, x2, . . . , xm and Y = y1, y2, . . . , yn and wish to find a maximum-length
common subsequence of X and Y. 

Algorithm Used: Dynamic Programing



About the Application:
The code requires few extra packages to be installed separately – numpy, sys, re 
Case1: For Global, Local and Semiglobal alignment

Input: Type of alignment, Strings for which sequences are to be aligned, Scoring values
Output: Alignment, Scoring matrix, and Score


Case2: of Longest Common Sub-Sequence

Input: Sequences for which LCS need to be found. 
Output: LCS, Length of the LCS and edit distance between the strings.


References:
0https://en.wikipedia.org/wiki/Sequence_alignment
1Algorithms on Strings, Trees, and Sequences - 11.2.1
2 https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
3Introduction to Computational Molecular biology- 3.2.2
4Algorithms on Strings, Trees,and Sequences - 3.2.3
5Introduction to Algorithms -Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, Clifford Stein- 15.4
