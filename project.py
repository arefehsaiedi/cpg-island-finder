"""
CpG Island Finder
MY CS50 Final Project
Simple tool to find CpG islands in DNA sequences from FASTA files
"""

import sys
import argparse
from fasta import read_fasta
from report import generate_report

def gc_percent(seq):
    """Return GC% of a DNA string"""
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) * 100

def cpg_ratio(seq):
    """Observed / expected CpG ratio in a DNA string"""
    n = len(seq)
    if n < 2:
        return 0.0

    c = seq.count("C")
    g = seq.count("G")

    # count CG dinucLeotides
    obs = 0
    for i in range(n - 1):
        if seq[i:i+2] == "CG":
            obs += 1

    if c == 0 or g == 0:
        return 0.0

    exp = (c * g) / n
    if exp == 0:
        return 0.0

    return obs / exp

def find_cpg_islands(seq, window_size=200, step=1, min_gc=50.0, min_ratio=0.6, do_merge=True):
        """
        Sliding window CpG island finder
        Returns a list of dicts with start, end, gc, ratio, lenght, and subsequence
        """
        n = len(seq)
        if n < window_size:
            return []

        window_size = int(window_size)
        step = max(1, int(step))

        raw_islands = []
        # fixed
        for start in range(0, n - window_size + 1, step):
               # check
            end = start + window_size
            win = seq[start:end]

            gc = gc_percent(win)   #change
            ratio = cpg_ratio(win)

            #debug print
            # print(f"window {start}-{end}: GC={gc:.1%}, Ratio={ratio:.1f}")

            if gc >= min_gc and ratio >= min_ratio:
                raw_islands.append((start, end, gc, ratio))
                #print(" -> Island found!")


        if do_merge and raw_islands:
            merged = merge_islands(raw_islands)
        else:
            merged = raw_islands

        result = []
        for st, en, gc, r in merged:
            islands_seq = seq[st:en]
            result.append({
                "start": st,
                "end": en,
                "length": en - st,
                "gc_content": gc,
                "cpg_ratio": r,
                "sequence": islands_seq
            })

        print(f"Total islands found: {len(result)}")
        return result

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



