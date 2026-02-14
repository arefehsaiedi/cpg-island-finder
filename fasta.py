"""
Simple FASTA reader for my CpG project
Just reads  a DNA sequence from a FASTA file
""" 

def read_fasta(filename):
    """Read DNA sequence from FASTA file"""
    seqs = {}
    name = None
    dna = []

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # save
                if name is not None:
                    seqs[name] = "".join(dna).upper()
                    # new
                name = line[1:]
                dna = []
            else:
                    # seq line
                    dna.append(line.replace(" ", ""))

    if name is not None:
        seqs[name] = "".join(dna).upper()

    return seqs  
   
        

def check_dna(dna):
    """Check if it's valid DNA (just ACGT)"""
    for base in dna:
        if base not in "ACGTN":
            return False
    return True
    
def dna_stats(dna):
    """Quick stats about the DNA sequence"""
    total = len(dna)
    a = dna.count("A")
    c = dna.count("C") 
    g = dna.count("G") 
    t = dna.count("T")     

    gc = (c + g) / total * 100 if total > 0 else 0

    return f"{total} bases, GC: {gc:.1f}%"
