from collections import defaultdict
from helper_math import calculate_results
from helper_math import perform_analysis
from helper_math import load_fasta_files

# SETTINGS
EVALUATE_RNA5S = True
EVALUATE_HOXA = True
EVALUATE_HOXB = True
EVALUATE_HOXC = True
EVALUATE_HOXD = True
DELETE_CACHE = False

# Evaluation code
def evaluate_gene(gene_name):
    fasta_files = load_fasta_files(f"data/{gene_name}/")
    if not fasta_files:
        print(f"No .fa files found in the directory data/{gene_name}/")
    else:
        calculate_results(fasta_files, gene_name)
        perform_analysis(gene_name)

# MAIN
if DELETE_CACHE:
    import shutil
    shutil.rmtree('cache')

if EVALUATE_RNA5S:
    evaluate_gene("RNA5S")

if EVALUATE_HOXA:
    evaluate_gene("HOXA")

if EVALUATE_HOXB:
    evaluate_gene("HOXB")

if EVALUATE_HOXC:
    evaluate_gene("HOXC")

if EVALUATE_HOXD:
    evaluate_gene("HOXD")


