#!/usr/bin/env python3
"""Diff two barcode-mapping runs.

  python compare_barcode.py OLD_MAPPING_DIR NEW_MAPPING_DIR
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

CHELON = r"tortoise|testudo|geochelone|astrochelys|stigmochelys|agrionemys|chelon"
FAMILY = ["ST10", "ST42", "ST43", "ST44", "ST10-complex"]


def load(d: Path):
    return pd.read_csv(d / "manifest_short_reads.csv")


def confident(df):
    return df[~df.confidence.isin(["unmapped", "ambiguous"])]


def chelonian_label(run: Path):
    """Which novel_N is the chelonian cluster in this run's reference set."""
    mf = pd.read_csv(run / "fasta" / "manifest.csv")
    nov = mf[mf.subtype.str.startswith("novel_", na=False)]
    for k, g in nov.groupby("subtype"):
        if len({h for h in g.host if re.search(CHELON, str(h), re.I)}) >= 3:
            return k
    return None


def row(lab, a, b, fmt="{:,}"):
    d = b - a
    print(f"  {lab:38s} {fmt.format(a):>10} -> {fmt.format(b):>10}   {d:+,}")


def main():
    old_d, new_d = Path(sys.argv[1]), Path(sys.argv[2])
    o, n = load(old_d), load(new_d)
    oc, nc = confident(o), confident(n)

    print(f"OLD {old_d}\nNEW {new_d}\n")
    print("BARCODE CORPUS")
    row("records retrieved (100-999 bp)", len(o), len(n))
    row("confident calls", len(oc), len(nc))
    print(f"  {'call rate':38s} {len(oc)/len(o)*100:>9.1f}% -> "
          f"{len(nc)/len(n)*100:>9.1f}%")
    for tier in ["high", "medium", "low", "ambiguous", "unmapped"]:
        row(f"  {tier}", int((o.confidence == tier).sum()),
            int((n.confidence == tier).sum()))

    print("\nST10 / ST42-44 FAMILY")
    for s in FAMILY:
        row(f"  {s}", int((o.assigned == s).sum()), int((n.assigned == s).sum()))
    if "separation_note" in n.columns:
        col = n.separation_note.fillna("")
        print(f"  {'collapsed to a complex':38s} {int(col.str.startswith('collapsed').sum()):>10,}")
        print(f"  {'unresolved -> ambiguous':38s} {int(col.str.startswith('unresolved').sum()):>10,}")

    print("\nNOVEL-CLUSTER CORROBORATION")
    ob = oc[oc.assigned.str.startswith("novel_", na=False)]
    nb = nc[nc.assigned.str.startswith("novel_", na=False)]
    row("  records landing in novel clusters", len(ob), len(nb))
    row("  clusters corroborated", ob.assigned.nunique(), nb.assigned.nunique())

    ochel = chelonian_label(Path(sys.argv[3])) if len(sys.argv) > 3 else "novel_3"
    nchel = chelonian_label(Path(sys.argv[4])) if len(sys.argv) > 4 else None
    if nchel:
        a = int((ob.assigned == ochel).sum())
        b = int((nb.assigned == nchel).sum())
        ah = int(((ob.assigned == ochel) & (ob.confidence == "high")).sum())
        bh = int(((nb.assigned == nchel) & (nb.confidence == "high")).sum())
        print(f"  chelonian ({ochel} -> {nchel}): {a} -> {b} records, "
              f"{ah} -> {bh} high-confidence")


if __name__ == "__main__":
    main()
