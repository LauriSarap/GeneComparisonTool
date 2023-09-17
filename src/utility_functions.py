from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import os
import pickle


def log_progress(progress, processed_records, total_records, elapsed_time):
    print(f"Progress: {progress:.2f}% completed | {processed_records}/{total_records} alignments | Time elapsed: {elapsed_time:.2f} minutes")


def write_aggregated_results_to_csv(new_file_name, df, maximum_possible_similarity):
    grouped_by_species_df = df.groupby('Species').agg({'Similarity': 'sum'}).reset_index()
    grouped_by_species_df['Similarity'] = (grouped_by_species_df['Similarity'] / maximum_possible_similarity) * 100
    grouped_by_species_df_sorted = grouped_by_species_df.sort_values(by='Similarity', ascending=False)
    grouped_by_species_df_sorted.to_csv(new_file_name, index=False)


def delete_csv_file(filename):
    if os.path.exists(filename):
        os.remove(filename)


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
    # print(best_alignment)
    #print(format_alignment(*alignments[0]))
    max_length = max(len(seq1), len(seq2))
    return (score / max_length) * 100, best_alignment

