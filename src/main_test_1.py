from Bio import SeqIO
from Bio.Align import PairwiseAligner
import matplotlib.pyplot as plt


file_path = "data/test_comparison_data.fa"


def read_fasta(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        # Extracting the species name from the record ID
        species_name = record.id.split("/")[0]
        sequences[species_name] = str(record.seq)
    return sequences


def compare_sequences_biopython(seq1, seq2):
    aligner = PairwiseAligner()
    score = aligner.score(seq1, seq2)
    length = len(seq1)
    similarity_percentage = (score / length) * 100
    return similarity_percentage


def main():

    # Reading the sequences from the FASTA file
    sequences = read_fasta(file_path)
    print(sequences.keys())

    # Printing the sequences
    # for species_name, sequence in sequences.items():
    #    print(f"{species_name}: {sequence}")
    #    print('\n')

    human_sequence = sequences['homo_sapiens']

    # Comparing other species against Homo sapiens
    similarity_with_human_biopython = {}
    for species, sequence in sequences.items():
        if species != 'homo_sapiens':
            similarity = compare_sequences_biopython(human_sequence, sequence)
            similarity_with_human_biopython[species] = similarity

    # Finding the species with the highest similarity to humans
    most_similar_species_biopython = max(similarity_with_human_biopython, key=similarity_with_human_biopython.get)
    most_similar_percentage_biopython = similarity_with_human_biopython[most_similar_species_biopython]

    print("Most similar species:", most_similar_species_biopython)
    print("Similarity percentage:", most_similar_percentage_biopython)

    # Plotting
    sorted_data = sorted(similarity_with_human_biopython.items(), key=lambda x: x[1])
    sorted_species = [item[0] for item in sorted_data]
    sorted_similarities = [item[1] for item in sorted_data]

    plt.figure(figsize=[10, 5])
    plt.bar(sorted_species, sorted_similarities, color=['blue', 'green', 'orange'])
    plt.xlabel('Species')
    plt.ylabel('Similarity Percentage (%)')
    plt.title('Genetic Similarity with Homo Sapiens (Human)')
    plt.ylim(min(sorted_similarities) - 5, max(sorted_similarities) + 5)
    plt.xticks(rotation=0)
    plt.grid(axis='y', linestyle='--')
    plt.show()


if __name__ == "__main__":
    main()

