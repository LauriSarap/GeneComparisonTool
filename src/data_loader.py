import os
import csv

existing_genes = set()
with open('evaluations.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    next(csvreader)
    for row in csvreader:
        gene_group, gene_name, _ = row
        existing_genes.add((gene_group, gene_name))

new_genes = set()
base_dir = 'data/'
for gene_group in os.listdir(base_dir):
    gene_group_path = os.path.join(base_dir, gene_group)
    if os.path.isdir(gene_group_path):
        for gene_name in os.listdir(gene_group_path):
            gene_name_path = os.path.join(gene_group_path, gene_name)
            if os.path.isdir(gene_name_path):
                new_genes.add((gene_group, gene_name))

genes_to_add = new_genes - existing_genes

with open('evaluations.csv', 'a', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    for gene_group, gene_name in genes_to_add:
        csvwriter.writerow([gene_group, gene_name, 'FALSE'])

print(f"Added {len(genes_to_add)} new genes to evaluations.csv.")
