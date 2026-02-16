"""
CpG Island Finder 
MY CS50 Final Project
Simple tool to find CpG islands in DNA sequences from FASTA files 
"""

import sys
import argparse 
from fasta import read_fasta
from analysis import find_cpg_islands
from report import generate_report

def main():
    parser = argparse.ArgumentParser(description="Find CpG islands in DNA!")
    parser.add_argument("input", help="FASTA file with DNA sequence")
    parser.add_argument("-w", "--window", type=int, default=200, help="window size (default 200)")
    parser.add_argument("-g", "--gc", type=float, default=0.5, help="minimum GC count (default 0.5)")
    parser.add_argument("-r", "--ratio", type=float, default=0.6, help="minimum CpG raito (default 0.6)")
    parser.add_argument("-o", "--output", help="save report to file")
    
    args = parser.parse_args()
    
    print("Loading DNA sequence...")

    sequences = read_fasta(args.input)

    if not sequences:
        print("NO sequences found!")
        return 1
    
    # take the first sequence for simplicity
    seq_name = list(sequences.keys()) [0]
    dna = sequences[seq_name]
    print(f'Loaded {len(dna)} bases from {seq_name}')

    print("Running analysis for CpG islands...")

    islands = find_cpg_islands(dna, args.window, step=1, min_gc=args.gc*100, min_ratio=args.ratio)  #this one need check
    print(f'Found {len(islands)} potential CpG island(s)!')

    report = generate_report(islands, seq_name, len(dna))

    if args.output:
        with open(args.output, "w") as f:
            f.write(report)
        print(f'Save to {args.output}')
    else:
        print("\n" + "="*40)
        print(report)
        print("="*40)

    return 0 

if __name__ == "__main__":
    sys.exit(main())

    

