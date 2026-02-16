"""
Simple CpG island finder functions
Used by project.py in my CS50 project
"""
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

def merge_islands(islands, max_gap=100):
    """
    Merge overlapping or close islands
    """
    if not islands:
        return []

    # sort by start
    islands.sort(key=lambda x: x[0])
    merged = []
    cur_start, cur_end, cur_gc, cur_ratio = islands[0]   # check

    for st, en, gc, r in islands[1:]:
        if st <= cur_end + max_gap:
            # extend island
            cur_end = en
            cur_gc = (cur_gc + gc) / 2
            cur_ratio = (cur_ratio + r) / 2
        else:
            merged.append((cur_start, cur_end, cur_gc, cur_ratio))  # check
            cur_start, cur_end, cur_gc, cur_ratio = st, en, gc, r      # check

    merged.append((cur_start, cur_end, cur_gc, cur_ratio))
    return merged

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

