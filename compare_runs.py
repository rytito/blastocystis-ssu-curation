#!/usr/bin/env python3
"""Diff two pipeline runs: acquisition, labels, separation statistics,
metadata composition and novel clusters.

  python compare_runs.py OLD_RUN_DIR NEW_RUN_DIR
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

UNK = {"unknown", "nan", "", "none", "fecal"}
CHELON = r"tortoise|testudo|geochelone|astrochelys|stigmochelys|agrionemys|chelon"


def load(run: Path):
    mf = pd.read_csv(run / "fasta" / "manifest.csv")
    ca = pd.read_csv(run / "clusters" / "cluster_assignments.csv")
    log = (run / "pipeline_log.txt").read_text(errors="replace")
    z = np.load(run / "distances" / "distance_matrix.npz", allow_pickle=True)
    return dict(run=run, mf=mf, ca=ca, log=log, D=z["D"], ids=z["ids"])


def logint(log, pat):
    m = re.search(pat, log)
    return int(m.group(1)) if m else None


def stats(r):
    mf, log = r["mf"], r["log"]
    D, ids = r["D"], r["ids"]
    lab = r["ca"].set_index("accession").subtype.reindex(ids).values
    is_st = np.array([bool(re.fullmatch(r"ST\d+", str(s))) for s in lab])
    sub = np.where(is_st)[0]
    L = lab[sub]
    iu = np.triu_indices(len(sub), 1)
    d = D[np.ix_(sub, sub)][iu]
    same = L[iu[0]] == L[iu[1]]
    intra, inter = d[same], d[~same]
    U = mannwhitneyu(intra, inter, alternative="two-sided").statistic
    delta = 2 * U / (len(intra) * len(inter)) - 1

    ctry = mf.country.astype(str).str.split(":").str[0].str.strip().replace(
        {"Czech Republic": "Czechia"})
    nov = mf[mf.subtype.str.startswith("novel_", na=False)]
    chel = None
    for k, g in nov.groupby("subtype"):
        if len({h for h in g.host if re.search(CHELON, str(h), re.I)}) >= 3:
            chel = (k, len(g))
    beaver = next((k for k, g in nov.groupby("subtype")
                   if any(re.search(r"castor\s+fiber", str(h), re.I) for h in g.host)), None)

    return dict(
        retrieved=logint(log, r"Fetching (\d+) records"),
        near_full=logint(log, r"Length filter.*?(\d+) records") or len(mf),
        dedup_removed=logint(log, r"Cross-ST deduplication: (?:\d+) . (?:\d+)")
        or logint(log, r"(\d+) cross-label duplicates"),
        final=len(mf),
        n_labels=mf.subtype.nunique(),
        n_st=sum(bool(re.fullmatch(r"ST\d+", str(s))) for s in mf.subtype.unique()),
        n_novel=nov.subtype.nunique(),
        n_clustered=len(nov),
        n_single=int((mf.subtype == "ST_unassigned").sum()),
        n_st_records=int(is_st.sum()),
        med_intra=float(np.median(intra)), med_inter=float(np.median(inter)),
        delta=float(delta), n_pairs=len(d),
        countries=len({c for c in ctry if c.lower() not in UNK}),
        africa=0,
        n_dated=int(mf.collection_date.astype(str).str.contains(r"(?:19|20)\d{2}").sum()),
        tech=mf.seq_tech.value_counts().to_dict(),
        st_counts=mf[mf.subtype.str.fullmatch(r"ST\d+", na=False)].subtype.value_counts().to_dict(),
        chel=chel, beaver=beaver,
        ctry_counts=ctry[ctry.str.lower() != "unknown"].value_counts().to_dict(),
    )


def row(label, a, b, fmt="{:,}"):
    if a is None or b is None:
        print(f"  {label:34s} {str(a):>12} -> {str(b):>12}")
        return
    delta = b - a
    sign = f"{delta:+,}" if isinstance(delta, int) else f"{delta:+.4f}"
    print(f"  {label:34s} {fmt.format(a):>12} -> {fmt.format(b):>12}   {sign}")


def main():
    old, new = (stats(load(Path(p))) for p in sys.argv[1:3])
    print(f"OLD {sys.argv[1]}\nNEW {sys.argv[2]}\n")

    print("ACQUISITION")
    for k, l in [("retrieved", "records retrieved"), ("near_full", "near-full-length"),
                 ("final", "final reference set")]:
        row(l, old[k], new[k])

    print("\nLABELS")
    for k, l in [("n_labels", "distinct labels"), ("n_st", "accepted/tentative subtypes"),
                 ("n_st_records", "records with an ST label"),
                 ("n_novel", "provisional novel clusters"),
                 ("n_clustered", "records in novel clusters"),
                 ("n_single", "ST_unassigned singletons")]:
        row(l, old[k], new[k])

    print("\nSEPARATION")
    row("median within-subtype", old["med_intra"], new["med_intra"], "{:.4f}")
    row("median between-subtype", old["med_inter"], new["med_inter"], "{:.4f}")
    row("Cliff's delta", old["delta"], new["delta"], "{:.4f}")
    row("pairwise comparisons", old["n_pairs"], new["n_pairs"])

    print("\nMETADATA")
    row("countries represented", old["countries"], new["countries"])
    row("records with a collection date", old["n_dated"], new["n_dated"])
    for t in sorted(set(old["tech"]) | set(new["tech"]),
                    key=lambda x: -new["tech"].get(x, 0)):
        row(f"  technology: {t}", old["tech"].get(t, 0), new["tech"].get(t, 0))

    print("\nFLAGGED CLUSTERS (identified by composition, not label)")
    print(f"  chelonian : {old['chel']} -> {new['chel']}")
    print(f"  Castor    : {old['beaver']} -> {new['beaver']}")

    print("\nSUBTYPES THAT CHANGED")
    keys = set(old["st_counts"]) | set(new["st_counts"])
    ch = [(k, old["st_counts"].get(k, 0), new["st_counts"].get(k, 0)) for k in keys]
    ch = sorted([c for c in ch if c[1] != c[2]], key=lambda c: -abs(c[2] - c[1]))
    if not ch:
        print("  none")
    for k, a, b in ch[:25]:
        print(f"  {k:8s} {a:5d} -> {b:5d}   {b - a:+d}")

    print("\nCOUNTRIES THAT CHANGED")
    keys = set(old["ctry_counts"]) | set(new["ctry_counts"])
    cc = [(k, old["ctry_counts"].get(k, 0), new["ctry_counts"].get(k, 0)) for k in keys]
    cc = sorted([c for c in cc if c[1] != c[2]], key=lambda c: -abs(c[2] - c[1]))
    if not cc:
        print("  none")
    for k, a, b in cc[:15]:
        print(f"  {k:22s} {a:5d} -> {b:5d}   {b - a:+d}")


if __name__ == "__main__":
    main()
