from Bio import SeqIO
from collections import defaultdict
from helper_math import calculate_similarity, calculate_similarity_with_print
import csv


# SETTINGS
UNALIGNED_DATA = True


results = defaultdict(dict)

fasta_files = []

if UNALIGNED_DATA:
    for i in range(1, 18):
        if i != 8:
            fasta_files.insert(i, f"data/Human_RNA5S{i}_orthologues.fa")

    for fasta_file in fasta_files:
        variant_name = fasta_file.split("_")[1]
        records = list(SeqIO.parse(fasta_file, "fasta"))

        homo_sapiens_seq = None
        for record in records:
            if "ENST" in record.id:
                homo_sapiens_seq = record.seq
                break

        for record in records:
            if "ENST" not in record.id:
                species = record.id
                similarity, best_alignment = calculate_similarity(homo_sapiens_seq, record.seq)
                results[variant_name][species] = {
                    'similarity': similarity,
                    'score': best_alignment[2],
                    'start': best_alignment[3],
                    'end': best_alignment[4]
                }

    print("Results by variant and species!")
    for variant_name, species_data in results.items():
        print(f"Variant: {variant_name}")
        for species, similarity in species_data.items():
            print(f"  Species: {species}, Similarity: {similarity}")

    with open('alignment_results.csv', 'w', newline='') as csvfile:
        fieldnames = ['Variant Name', 'Species', 'Similarity', 'Score', 'Start', 'End']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for variant_name, species_data in results.items():
            for species, data in species_data.items():
                writer.writerow({
                    'Variant Name': variant_name,
                    'Species': species,
                    'Similarity': data['similarity'],
                    'Score': data['score'],
                    'Start': data['start'],
                    'End': data['end']
                })

