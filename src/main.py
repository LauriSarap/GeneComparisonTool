from collections import defaultdict
from helper_math import calculate_results
from helper_math import perform_analysis
from helper_math import load_fasta_files
import csv

# SETTINGS

evaluation_settings = {}
settings = {}

with open('evaluations.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    next(csvreader)
    for row in csvreader:
        category, gene, evaluate = row
        evaluation_settings[gene] = (evaluate.lower() == 'true', category)

with open('settings.txt', 'r') as settings_file:
    for line in settings_file:
        if "=" in line:
            key, value = line.strip().split('=', 1)
            key = key.strip()
            value = value.strip()
            if value.lower() == 'true':
                settings[key] = True
            elif value.lower() == 'false':
                settings[key] = False
            elif value.isdigit():
                settings[key] = int(value)
            else:
                settings[key] = value

for key, value in settings.items():
    print(f"{key}: {value}")

LOGGING_INTERVAL = settings.get('LOGGING_INTERVAL', 10)


# Evaluation code
def evaluate_gene(gene_name, gene_path):
    fasta_files = load_fasta_files(gene_path)
    if not fasta_files:
        print(f"No .fa files found in the directory data/{gene_name}/")
    else:
        calculate_results(fasta_files, gene_name, gene_path, LOGGING_INTERVAL)
        perform_analysis(gene_name, gene_path)


for gene, (should_evaluate, category) in evaluation_settings.items():
    if should_evaluate:
        gene_path = f'data/{category}/{gene}/'
        evaluate_gene(gene, gene_path)


