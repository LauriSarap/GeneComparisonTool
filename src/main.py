from collections import defaultdict
from helper_math import calculate_results
from helper_math import perform_analysis

# SETTINGS
EVALUATE_RNA5S = False
EVALUATE_HOXA = True

# MAIN
if EVALUATE_RNA5S:
    RNA5S_fasta_files = []
    for i in range(1, 18):
        if i != 8:
            RNA5S_fasta_files.insert(i, f"data/RNA5S/Human_RNA5S{i}_orthologues.fa")

    calculate_results(RNA5S_fasta_files, "RNA5S")
    perform_analysis("RNA5S")

if EVALUATE_HOXA:
    HOXA_fasta_files = []
    for i in range(1, 14):
        if i != 8 and i != 11:
            HOXA_fasta_files.insert(i, f"data/HOXA/Human_HOXA{i}_orthologues.fa")

    calculate_results(HOXA_fasta_files, "HOXA")
    perform_analysis("HOXA")

