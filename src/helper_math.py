from Bio import pairwise2, SeqIO
from Bio.pairwise2 import format_alignment
from collections import defaultdict
import csv
import pandas as pd
import os
import re
import pickle
import time

ignore_list = ['MGP_CBAJ_', 'MGP_PWKPhJ_', 'MGP_C57BL6NJ_',
               'MGP_AKRJ_', 'MGP_DBA2J_', 'MGP_NZOHlLtJ_',
               'MGP_WSBEiJ_', 'MGP_LPJ_', 'MGP_CASTEiJ_',
               'MGP_BALBcJ_', 'MGP_129S1SvImJ_', 'MGP_FVBNJ_',
               'MGP_C3HHeJ_', 'MGP_NODShiLtJ_']


def load_fasta_files(directory):
    fasta_files = []
    for filename in os.listdir(directory):
        if filename.endswith(".fa"):
            fasta_files.append(os.path.join(directory, filename))
    return fasta_files


def calculate_similarity(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]
    score = best_alignment[2]
    #print(best_alignment)
    #print(format_alignment(*alignments[0]))
    max_length = max(len(seq1), len(seq2))
    return (score / max_length) * 100, best_alignment


def calculate_results(fasta_files, gene_name, gene_path, logging_interval=10, minimum_similarity=90):

    results = defaultdict(dict)
    if not os.path.exists('cache'):
        print(f"Creating directory cache!")
        os.makedirs('cache')

    cache_dir = f'cache/{gene_name}/'

    # Perform alignment
    print(f'Calculating results for {gene_name}...')
    for fasta_file in fasta_files:
        start_time = time.time()
        variant_name = fasta_file.split("_")[1]
        cache_file = os.path.join(cache_dir, f"{variant_name}.pkl")

        # Create cache directory if not exists
        if not os.path.exists(cache_dir):
            print(f"Creating directory {cache_dir}")
            os.makedirs(cache_dir)

        # Load from alignment from cache if exists
        if os.path.exists(cache_file):
            print(f'Loading cached alignment for {variant_name}..')
            with open(cache_file, 'rb') as f:
                results[variant_name] = pickle.load(f)
            continue

        print("Processing file: " + variant_name + "..")
        records = list(SeqIO.parse(fasta_file, "fasta"))

        homo_sapiens_seq = None
        for record in records:
            if re.match(r'^ENS[T|P]\d+$', record.id):  # Homo sapiens
                #print("Found homo sapiens for: " + variant_name)
                #print(record.id)
                homo_sapiens_seq = record.seq
                break

        if homo_sapiens_seq is not None:
            total_records = sum(1 for record in records if not re.match(r'^ENS[T|P]\d+$', record.id))
            processed_records = 0

            print(f"Progress: 0% completed | {processed_records}/{total_records} alignments | Time elapsed: 0 minutes")

            for record in records:
                if not re.match(r'^ENS[T|P]\d+$', record.id):  # Exclude Homo sapiens
                    processed_records += 1

                    if any(record.id.startswith(prefix) for prefix in ignore_list):
                        print(f"Ignoring {record.id}")
                        continue

                    species = record.id
                    similarity, best_alignment = calculate_similarity(homo_sapiens_seq, record.seq)
                    results[variant_name][species] = {
                        'similarity': similarity,
                        'score': best_alignment[2],
                        'start': best_alignment[3],
                        'end': best_alignment[4]
                    }

                    progress = (processed_records / total_records) * 100
                    if processed_records % logging_interval == 0:
                        elapsed_time = (time.time() - start_time) / 60
                        print(f"Progress: {progress:.2f}% completed | {processed_records}/{total_records} alignments | Time elapsed: {elapsed_time:.2f} minutes")
        else:
            print(f"No Homo sapiens sequence found in {fasta_file}. Skipping.")

        # Cache the alignment results
        print(f"Caching alignment for {variant_name}...")
        with open(cache_file, "wb") as f:
            pickle.dump(results[variant_name], f)

        total_time_elapsed = (time.time() - start_time) / 60
        print(f'Finished processing {fasta_file} | Total time: {total_time_elapsed:.2f} minutes')

    delete_csv_file(f'{gene_path}{gene_name}_alignment_results.csv')

    # Map animal names to species prefixes
    species_mapping = {}
    with open("data/species_prefixes.txt", "r") as f:
        next(f)  # Skip the header
        for line in f:
            prefix, name = line.strip().split("\t")
            species_mapping[prefix] = name

    # Save alignment scores to .csv
    print(f'Saving {gene_name} results to .csv..')
    with open(f"{gene_path}{gene_name}_alignment_results.csv", 'w', newline='') as csvfile:
        fieldnames = ['Variant Name', 'Species', 'Similarity', 'Score', 'Start', 'End']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for variant_name, species_data in results.items():
            for species, data in species_data.items():
                species_prefix = max((k for k in species_mapping.keys() if species.startswith(k)), key=len, default=None)
                species_name = species_mapping.get(species_prefix, species)
                if data['similarity'] >= minimum_similarity:
                    writer.writerow({
                        'Variant Name': variant_name,
                        'Species': species_name,
                        'Similarity': data['similarity'],
                        'Score': data['score'],
                        'Start': data['start'],
                        'End': data['end']
                    })
    print(f'Finished saving {gene_name} results to .csv!')


def perform_analysis(gene_name, gene_path):
    delete_csv_file(f"{gene_path}{gene_name}_aggregated_similarity.csv")
    print(f'Performing analysis for {gene_name}..')

    df = pd.read_csv(f"{gene_path}{gene_name}_alignment_results.csv")

    # Filter out duplicate species for each variant based on the highest similarity
    df = df.loc[df.groupby(['Variant Name', 'Species'])['Similarity'].idxmax()]

    maximum_possible_similarity = len(df['Variant Name'].unique()) * 100

    grouped_by_species_df = df.groupby('Species').agg({'Similarity': 'sum'}).reset_index()
    grouped_by_species_df['Similarity'] = (grouped_by_species_df['Similarity'] / maximum_possible_similarity) * 100
    grouped_by_species_df_sorted = grouped_by_species_df.sort_values(by='Similarity', ascending=False)
    grouped_by_species_df_sorted.to_csv(f"{gene_path}{gene_name}_aggregated_similarity.csv", index=False)

    print(f'Finished performing analysis for {gene_name}!')


def perform_analysis_for_group(gene_group_name, gene_group_path):
    delete_csv_file(f"{gene_group_path}{gene_group_name}_aggregated_similarity.csv")
    print(f'Performing analysis for gene group {gene_group_name}..')

    genes_to_evaluate = []
    with open('evaluations.csv', 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader)
        for row in csvreader:
            group, gene, evaluate = row
            if group == gene_group_name:
                genes_to_evaluate.append(gene)

    maximum_possible_similarity = len(genes_to_evaluate) * 100

    aggregated_df = pd.DataFrame()

    for gene in genes_to_evaluate:
        gene_csv_path = os.path.join(gene_group_path, gene, f"{gene}_aggregated_similarity.csv")
        gene_df = pd.read_csv(gene_csv_path)
        aggregated_df = pd.concat([aggregated_df, gene_df])

    aggregated_analysis_df = aggregated_df.groupby('Species').agg({'Similarity': 'sum'}).reset_index()
    aggregated_analysis_df['Similarity'] = (aggregated_analysis_df['Similarity'] / maximum_possible_similarity) * 100
    aggregated_analysis_df_sorted = aggregated_analysis_df.sort_values(by='Similarity', ascending=False)
    aggregated_analysis_df_sorted.to_csv(f"{gene_group_path}{gene_group_name}_group_aggregated_similarity.csv", index=False)


def perform_analysis_for_all_groups(gene_groups, basedir='data/'):
    print(f'Performing analysis for all gene groups..')
    overall_aggregated_df = pd.DataFrame()

    for gene_group in gene_groups:
        group_aggregated_df = pd.read_csv(f"{basedir}{gene_group}/{gene_group}_group_aggregated_similarity.csv")
        overall_aggregated_df = pd.concat([overall_aggregated_df, group_aggregated_df])

    maximum_possible_similarity = len(gene_groups) * 100

    overall_aggregated_analysis_df = overall_aggregated_df.groupby('Species').agg({'Similarity': 'sum'}).reset_index()
    overall_aggregated_analysis_df['Similarity'] = (overall_aggregated_analysis_df['Similarity'] / maximum_possible_similarity) * 100
    overall_aggregated_csv_path = os.path.join(basedir, 'overall_aggregated_analysis.csv')
    overall_aggregated_analysis_df.to_csv(overall_aggregated_csv_path, index=False)


def delete_all_csv_files(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            os.remove(os.path.join(directory, filename))


def delete_csv_file(filename):
    if os.path.exists(filename):
        os.remove(filename)


