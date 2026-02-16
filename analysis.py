"""
Simple CpG island finder functions
Used by project.py in my CS50 project
"""

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



