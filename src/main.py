from collections import defaultdict
from helper_math import calculate_results
from helper_math import perform_analysis
from helper_math import load_fasta_files
from helper_math import perform_analysis_for_group
from helper_math import perform_analysis_for_all_groups
import csv

# SETTINGS

evaluation_settings = {}
group_evaluation_settings = {}
settings = {}

with open('evaluations.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    next(csvreader)
    for row in csvreader:
        category, gene, evaluate = row
        evaluation_settings[gene] = (evaluate.lower() == 'true', category)

with open('group_evaluations.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    next(csvreader)
    for row in csvreader:
        category, evaluate = row
        group_evaluation_settings[category] = (evaluate.lower() == 'true')

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

LOGGING_INTERVAL = settings.get('LOGGING_INTERVAL', 10)
MINIMUM_SIMILARITY_PERCENTAGE = settings.get('MINIMUM_SIMILARITY_PERCENTAGE', 90)
PERFORM_GENE_EVALUATIONS = settings.get('PERFORM_GENE_EVALUATIONS', True)
PERFORM_GENE_GROUP_EVALUATIONS = settings.get('PERFORM_GENE_GROUP_EVALUATIONS', True)
PERFORM_OVERALL_EVALUATIONS_BY_GENE_VARIANT = settings.get('PERFORM_OVERALL_EVALUATIONS_BY_GENE_VARIANT', True)
PERFORM_OVERALL_EVALUATIONS_BY_GENE = settings.get('PERFORM_OVERALL_EVALUATIONS_BY_GENE', True)
PERFORM_OVERALL_EVALUATIONS_BY_GENE_GROUP = settings.get('PERFORM_OVERALL_EVALUATIONS_BY_GENE_GROUP', True)


# Evaluation code
def evaluate_gene(gene_name, gene_path):
    fasta_files = load_fasta_files(gene_path)
    if not fasta_files:
        print(f"No .fa files found in the directory data/{gene_name}/")
    else:
        calculate_results(fasta_files, gene_name, gene_path, LOGGING_INTERVAL, MINIMUM_SIMILARITY_PERCENTAGE)
        perform_analysis(gene_name, gene_path)


if PERFORM_GENE_EVALUATIONS:
    for gene, (should_evaluate, category) in evaluation_settings.items():
        if should_evaluate:
            gene_path = f'data/{category}/{gene}/'
            evaluate_gene(gene, gene_path)

if PERFORM_GENE_GROUP_EVALUATIONS:
    for group, should_evaluate in group_evaluation_settings.items():
        if should_evaluate:
            group_path = f'data/{group}/'
            perform_analysis_for_group(group, group_path)

if PERFORM_OVERALL_EVALUATIONS_BY_GENE_VARIANT:
    gene_variants = []
    for gene, (should_evaluate, category) in evaluation_settings.items():
        if should_evaluate:
            gene_variants.append(gene)

if PERFORM_OVERALL_EVALUATIONS_BY_GENE_GROUP:
    gene_groups = []
    for group, evaluations in group_evaluation_settings.items():
        gene_groups.append(group)
        print(group)

    perform_analysis_for_all_groups(gene_groups)

