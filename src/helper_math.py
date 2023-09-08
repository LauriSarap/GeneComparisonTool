from Bio import pairwise2
from Bio.pairwise2 import format_alignment


def calculate_similarity(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]
    score = best_alignment[2]
    #print(best_alignment)
    #print(format_alignment(*alignments[0]))
    max_length = max(len(seq1), len(seq2))
    return (score / max_length) * 100, best_alignment


def calculate_similarity_with_print(seq1, seq2):
    # Assuming both sequences are of the same length

    print("Homo sapiens sequence:")
    print(seq1)
    print("Other sequence:")
    print(seq2)

    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (matches / len(seq1)) * 100
