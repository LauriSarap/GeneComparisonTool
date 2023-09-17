# Homeobox Gene Comparison and Conservation Tool

## Overview

This repository contains the Python code and associated files used for analyzing the conservation of homeobox genes across various non-primate species compared to Homo sapiens. 
The research aims to answer the question: **"How similar are homeobox genes in various non-primates to those in Homo sapiens, and which species show the most significant levels of conservation?"**

## Features

- Sequence alignment done via Biopython pairwise aligner
- Percent identity calculation
- Data extraction from the Ensembl database
- CSV output containing percent identity scores for each gene
- CSV output for each Homeobox gene class
- Sequences prealigned for Homeobox genes coding sequences (Pseudogenes not included)

## Requirements

- Python 3.9
- Biopython
- pandas
- NumPy

## Usage

1. Using console clone the repository
```
git clone https://github.com/LauriSarap/GeneComparisonTool.git
```


2. Navigate to the repository folder
```
cd GeneComparisonTool
```
  

3. Install the required Python packages
```
pip install biopython
pip install pandas
```
  

4. Run the main python script
```
cd src
python main.py
```


## Files

- `main.py`: The main Python script for running the analysis.
- `src/data/`: Folder containing the raw FASTA files and where .csv files are created
- `src/cache/`: Folder where alignments are cached

## Limitations

- The code is tailored to work with pre-computed potential orthologous genes from the Ensembl database.
- The current scoring system for sequence alignment is basic, giving a score of +1 for a match and 0 for a mismatch.

## Acknowledgments

This research was conducted using data from reputable sources, such as the HUGO Gene Nomenclature Committee (HGNC) and Ensembl databases.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
