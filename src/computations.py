from Bio import SeqIO
from collections import defaultdict
from utility_functions import *
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


def calculate_alignment_results(fasta_files, gene_name, gene_path, logging_interval=10, minimum_similarity=90):

    results = defaultdict(dict)

    # Cache directory
    if not os.path.exists('cache'):
        print(f"Creating directory cache!")
        os.makedirs('cache')
    cache_dir = f'cache/{gene_name}/'

    # Perform alignment
    for fasta_file in fasta_files:
        start_time = time.time()
        variant_name = fasta_file.split("_")[1]
        cache_file = os.path.join(cache_dir, f"{variant_name}.pkl")

        # Gene cache directory
        if not os.path.exists(cache_dir):
            print(f"Creating directory {cache_dir}")
            os.makedirs(cache_dir)

        # Load  alignment from cache if exists
        if os.path.exists(cache_file):
            print(f'Loading cached alignment for {variant_name}..')
            with open(cache_file, 'rb') as f:
                results[variant_name] = pickle.load(f)
            continue

        print("Computing alignment for: " + variant_name + "..")
        records = list(SeqIO.parse(fasta_file, "fasta"))

        homo_sapiens_seq = None
        for record in records:
            if re.match(r'^ENS[T|P]\d+$', record.id):  # Homo sapiens
                homo_sapiens_seq = record.seq
                break

        if homo_sapiens_seq is not None:
            total_records = sum(1 for record in records if not re.match(r'^ENS[T|P]\d+$', record.id))
            processed_records = 0
            log_progress(0, processed_records, total_records, 0)

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
                        log_progress(progress, processed_records, total_records, elapsed_time)
        else:
            print(f"No Homo sapiens sequence found in {fasta_file}. Skipping.")

        # Cache the alignment results
        print(f"Caching alignments for {variant_name}...")
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


def calculate_highest_similarity(gene_name, gene_path):
    delete_csv_file(f"{gene_path}{gene_name}_aggregated_similarity.csv")
    print(f'Calculating highest similarity for {gene_name}..')

    df = pd.read_csv(f"{gene_path}{gene_name}_alignment_results.csv")

    # Filter out duplicate species for each variant based on the highest similarity
    df = df.loc[df.groupby(['Variant Name', 'Species'])['Similarity'].idxmax()]

    maximum_possible_similarity = len(df['Variant Name'].unique()) * 100
    write_aggregated_results_to_csv(f"{gene_path}{gene_name}_aggregated_similarity.csv", df, maximum_possible_similarity)


def calculate_highest_similarity_gene_group(gene_group_name, gene_group_path):
    delete_csv_file(f"{gene_group_path}{gene_group_name}_aggregated_similarity.csv")
    print(f'Calculating highest similarity for group {gene_group_name}..')
    df = pd.DataFrame()

    genes_to_evaluate = []
    with open('evaluations.csv', 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader)
        for row in csvreader:
            group, gene, evaluate = row
            if group == gene_group_name:
                genes_to_evaluate.append(gene)
                
    for gene in genes_to_evaluate:
        gene_csv_path = os.path.join(gene_group_path, gene, f"{gene}_aggregated_similarity.csv")
        gene_df = pd.read_csv(gene_csv_path)
        df = pd.concat([df, gene_df])
        
    maximum_possible_similarity = len(genes_to_evaluate) * 100
    write_aggregated_results_to_csv(f"{gene_group_path}{gene_group_name}_group_aggregated_similarity.csv", df, maximum_possible_similarity)


def calculate_highest_similarity_overall_based_on_groups(gene_groups, basedir='data/'):
    print(f'Calculating highest similarity based on groups..')
    df = pd.DataFrame()
    
    for gene_group in gene_groups:
        group_df = pd.read_csv(f"{basedir}{gene_group}/{gene_group}_group_aggregated_similarity.csv")
        df = pd.concat([df, group_df])

    maximum_possible_similarity = len(gene_groups) * 100
    write_aggregated_results_to_csv(f"{basedir}overall_aggregated_similarity.csv", df, maximum_possible_similarity)


