CpG Island Finder - CS50 Final Project

This repo holds my final project for CS50:
Introduction to Programming with Python. 
It's a simple Python tool for finding CpG islands in DNA sequences from FASTA files.
The program reads a DNA sequence, uses a sliding window to calculate GC content and CpG ratios, spots potential islands, merges close ones, and creates a clear text report.


Project Overview
What it does:
It scans FASTA files for CpG islands with a sliding window. It finds regions above GC and CpG ratio thresholds, then merges nearby islands to avoid splits. 

Key features:
- Command-line options for window size (default: 200 bp), GC thresholds (0.5), CpG ratio (0.6), and output file.
- Handles test sequences and real data like BRCA1.
- Built with four main Python files for modularity.
- Produces easy-to-read reports with stats, island tables, and sequence snippets.
- Uses only Python's standard library (no extra install needed).

Project Structure
cpg-island-finder/
- main.py      # Command-line entry point
- fasta.py     # FASTA parsing
- analysis.py  # CpG detection logic
- report.py    # Report creation
- README.md    # This doc
- requirements.txt  # No dependencies
- usage.txt    # Quick examples
- test_fasta   # Test files
- example_outputs   # Sample reports

How the Algorithm Works
- Slide a fixed window along the DNA sequence.
- Calculate GC% and observed/expected CpG ratio per window.
- Filter by user thresholds.
- Merge overlapping or close islands (max gap: 100 bp).
- Output islands with positions, length, GC%, ratio, and preview.

Quick Demo
# Test sequence (CpG-rich) - finds one 490 bp island
python main.py test_fasta/test_cpg_rich.fasta -w 50 -g 0.5 -r 0.6

# BRCA1 gene (real NCBI data) - find multiple islands
python main.py test_fasta/brca1_sequence.fasta -w 50 -g 0.5 -r 0.6 -o brca1_report.txt

Test Output
===================================
===================================
CpG islands finder
===================================
sequence : NM_007294.4 Homo sapiens BRCA1 DNA repair associated (BRCA1), transcript variant 1, mRNA
length : 7088 bp
time : 2026-02-10 18:24:59:

summary
------
number of islands : 3
total islands bp    : 197 bp
coverage            : 2.78%
avg length          : 65.7 bp
avg GC              : 50.00%
avg CpG ratio       : 0.652

islands (positions are 0-based):
--------------------------------
id  start   end    len  GC%   ratio
--  -----   -----  ---- ----- -----
1   508     558     50     50.0    0.667
2   3821    3890    69     50.0    0.648
3   4707    4785    78     50.0    0.641

Setup and Usage 
Just need Python (no install required).

pip install -r requirements.txt   # No external packages required
python main.py input.fasta [-w 300] [-g 0.55] [-r 0.7] [-o output.txt]

Edge cases handled
- Empty or bad FASTA files.
- Invalid DNA letters.
- No CpG islands found.
- Overlapping or close regions.

What I Learned
- Basic Bioinformatics: CpG islands, GC content, ratios.
- Sliding window techniques for sequences.
- Modular Python code that's easy to maintain.
- Parsing FASTA and checking bio data.
- Creating clear scientific reports.
- Testing with fake and real sequences.

Future Ideas
- Add plots for GC and islands (using matplotlib).
- Multi-FASTA support.
- JSON/CSV outputs.
- Speed up for big genomes.

Acknowledgments
- CS50 (Harvard) for the programming basics.
- Gardiner-Garden & Frommer for CpG criteria.
- NCBI for free genomic data.
- Python community for docs and help.

This project let me use CS50 skills on a real bioinformatics task. Building the detector from scratch helped with modularity, CLI tools, and algorithms, plus an introdution to genomics.
It shows how Python makes simple, reliable tools for bio data analysis in a class setting.