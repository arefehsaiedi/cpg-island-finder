![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)

![Status](https://img.shields.io/badge/status-completed-brightgreen)

![License](https://img.shields.io/badge/license-MIT-lightgrey)



# CpG Island Finder 

### A Python Tool for Introductory Genomic Sequence Analysis 



CpG Island Finder is a Python-based bioinformatics tool developed as the final project for Harvard University’s CS50: Introduction to Programming with Python. 


The program analyses DNA sequences in FASTA format and identifies potential CpG islands using a sliding-window approach based on GC content and observed/expected CpG ratios. The project combines programming fundamentals with concepts from genomics, epigenetics, and molecular genetics.


# Biological Background 

CpG islands are genomic regions enriched in CpG dinucleotides and characterized by elevated GC content. They are commonly associated with gene promoters and play important roles in gene regulation, DNA methylation, cancer biology, and epigenetics.


This project was designed as an educational introduction to computational genomics, demonstrating how Python can be applied to biological sequence analysis.



# Repository Highlights 

- Developed as the final project for Harvard University’s CS50: Introduction to Programming with Python. 

- Implements a complete CpG island detection workflow from FASTA parsing to report generation. 

- Tested on both synthetic datasets and real human genomic sequences (BRCA1).

- Built entirely with Python standard libraries.

- Developed and documented on GitHub using version control practices.



# Project Overview

What it does:

It scans FASTA files for CpG islands with a sliding window. It finds regions above GC and CpG ratio thresholds, then merges nearby islands to avoid splits.



# Key features

- Command-line options for window size (default: 200 bp), GC thresholds (0.5), CpG ratio (0.6), and output file.

- Handles test sequences and real data like BRCA1.

- Built with four main Python files for modularity.

- Produces easy-to-read reports with stats, island tables, and sequence snippets.

- Uses only Python's standard library (no extra install needed).



# Project Structure

cpg-island-finder/

- project.py   # Command-line entry point and core CpG analysis function

- fasta.py     # FASTA parsing

- analysis.py  # Helper functions (optional)

- report.py    # Report creation

- README.md    # This doc

- requirements.txt  # No dependencies

- usage.txt    # Quick examples

- test_fasta   # Test files

- example_outputs   # Sample reports

- test_project.py  # Simple unit tests for my function



# How the Algorithm Works

- Slide a fixed window along the DNA sequence.

- Calculate GC% and observed/expected CpG ratio per window.

- Filter by user thresholds.

- Merge overlapping or close islands (max gap: 100 bp).

- Output islands with positions, length, GC%, ratio, and preview.



# Quick Demo

# Test sequence (CpG-rich) - finds one 490 bp island

python project.py test_fasta/test_cpg_rich.fasta -w 50 -g 0.5 -r 0.6


# BRCA1 gene (real NCBI data) - find multiple islands

python project.py test_fasta/brca1_sequence.fasta -w 50 -g 0.5 -r 0.6 -o brca1_report.txt


<img width="800" alt="First few lines of BRCA1 FASTA file showing sequence header and DNA bases" src="https://github.com/user-attachments/assets/f51f5e26-80fa-486b-8ff0-c4729cd17d3e" />

# Test Output

<img width="800" alt="brca1_output" src="https://github.com/user-attachments/assets/40fe6d82-6f01-4c6c-8198-3915406db097" />


# Setup and Usage

Just need Python (no install required).

python project.py input.fasta [-w 300] [-g 0.55] [-r 0.7] [-o output.txt]



# Edge cases handled

- Empty or bad FASTA files.

- Invalid DNA letters.

- No CpG islands found.

- Overlapping or close regions.



# What I Learned

- Applying Python programming to biological sequence analysis.

- Working with genomic data in FASTA format.

- Implementing sliding-window algorithms for feature detection.

- Designing modular and maintainable scientific software.

- Generating reproducible analysis reports.

- Connecting concepts from genomics, epigenetics, and molecular genetics with computational methods.



# Future Ideas

- Add plots for GC and islands (using matplotlib).

- Multi-FASTA support.

- JSON/CSV outputs.

- Speed up for big genomes.



# Acknowledgments

- CS50 team (Harvard) - for the programming basics.

- Gardiner-Garden & Frommer - for the classic CpG island criteria.

- NCBI - for providing free genomic data (BRCA1 sequence).

- Python community - for excellent documentation and support.


This project allowed me to combine my interest in genetics with programming while gaining hands-on experience in bioinformatics and genomic sequence analysis. It represents my first step toward applying computational approaches to molecular biology research.
