"""
Simple tests for CpG Island Finder
Just run this file to test everything
"""
from project import gc_percent, cpg_ratio, find_cpg_islands

def test_gc_percent():
    #empty seq
    assert gc_percent("") == 0.0

    # only A and T
    assert gc_percent("ATATAT") == 0.0

    # Mixed seq
    assert gc_percent("GATC") == 50.0
    assert gc_percent("GGCC") == 100.0
    assert gc_percent("AATT") == 0.0

    print("gc_percent tests passed")

def test_cpg_ratio():
    # too short
    assert cpg_ratio("") == 0.0
    assert cpg_ratio("A") == 0.0

    # no C or no G
    assert cpg_ratio("ATATAT") == 0.0

    # C and G present but no "CG" dinucleotide
    result =  cpg_ratio("ACGT")
    assert result < 0.001

    # CpG_rich seq
    assert cpg_ratio("CGCG") > 0
    assert cpg_ratio("ATCGATCG") > 0

    print("cpg_ratio tests passed")

def test_find_cpg_islands():
    # empty seq
    assert find_cpg_islands("") == []

    # shorter than window size
    short = "ATGC" * 10 #40 bp
    assert find_cpg_islands(short, window_size=50) == []

    # CpG_rich seq
    rich = "CG" * 100 #200 bp
    islands1 = find_cpg_islands(rich, window_size=100, min_gc=50.0, min_ratio=0.6)
    assert len(islands1) > 0

    # AT-rich seq
    poor = "AT" * 200 #400 bp
    islands2 = find_cpg_islands(poor, window_size=100, min_gc=50.0, min_ratio=0.6)
    assert len(islands2) == 0

    # check structure of 1 island
    if islands1:
        isl = islands1[0]
        assert "start" in isl
        assert "end" in isl
        assert "length" in isl
        assert "gc_content" in isl
        assert "cpg_ratio" in isl
        assert "sequence" in isl

    print("find_cpg_islands tests passed")

def main():
    print("Running tests for CpG Island Finder...\n")

    test_gc_percent()
    test_cpg_ratio()
    test_find_cpg_islands()

    print("\nAll tests passed.")

if __name__ == "__main__":
    main()
