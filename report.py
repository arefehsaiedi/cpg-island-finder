"""
The simple text report for CpG islands
Used by main.py 
"""

from datetime import datetime

def generate_report(islands, seq_name, seq_len):
    """
    islands: list of dicts from find_cpg_islands
    seq_name: name/id of sequence
    seq_len: length of sequence
    """
    lines = []

    # header (nothing fancy)
    lines.append("=" * 35)
    lines.append("CpG islands finder")
    lines.append("=" * 35)
    lines.append(f'sequence : {seq_name}')
    lines.append(f'length : {seq_len} bp')
    lines.append(f'time : {datetime.now().strftime('%Y-%m-%d %H:%M:%S:')}')
    lines.append("")

    # how many islands
    lines.append("summary")
    lines.append("------")
    lines.append(f'number of islands : {len(islands)}')

    if islands:
        total_bp = 0
        gc_sum = 0.0
        ratio_sum = 0.0

        for isl in islands:
            total_bp += isl["length"]
            gc_sum += isl["gc_content"]
            ratio_sum += isl["cpg_ratio"]

        cov = (total_bp / seq_len) * 100 if seq_len > 0 else 0
        avg_len = total_bp / len(islands)
        avg_gc = gc_sum / len(islands)
        avg_ratio = ratio_sum / len(islands)

        lines.append(f'total islands bp    : {total_bp} bp')
        lines.append(f'coverage            : {cov:.2f}%')
        lines.append(f'avg length          : {avg_len:.1f} bp')
        lines.append(f'avg GC              : {avg_gc:.2f}%')
        lines.append(f'avg CpG ratio       : {avg_ratio:.3f}')
    lines.append("")

    # details (kind of rough table)
    if islands:
        lines.append("islands (positions are 0-based):")
        lines.append("--------------------------------")
        lines.append("id  start   end    len  GC%   ratio")
        lines.append("--  -----   -----  ---- ----- -----")

        for i, isl in enumerate(islands, start=1):
            s = isl["start"]
            e = isl["end"]
            l = isl["length"]
            gc = isl["gc_content"]
            r = isl["cpg_ratio"]
            lines.append(
                f"{i:<2}  {s:<6}  {e:<6}  {l:<5}  {gc:<6.1f}  {r:<.3f}"
            )

        lines.append("")
        lines.append("island sequence (first 50 bp):")
        lines.append("------------------------------")
        for i, isl in enumerate(islands, start=1):
            subseq = isl["sequence"]
            if len(subseq) > 50:
                lines.append(f"{i}. {subseq}")

    # if no islands just leave summary part
    if not islands:
        lines.append("no CpG islands found :(")

    return "\n".join(lines)       

