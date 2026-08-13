"""
Map short (barcode-region) Blastocystis SSU rRNA records to the curated
reference set
============================================================================

The main pipeline (blastocystis_rRNA_pipeline_v5.py) deliberately excludes
the ~16,000 records shorter than 1000 bp — mostly the ~600 bp Scicluna
barcode region that dominates surveillance sequencing. Those records cannot
support novel-ST discovery, but they CAN be subtype-called against a curated
full-length reference set. This script does exactly that:

  1. Re-issues the same broad NCBI query as the main pipeline.
  2. Keeps the records the main pipeline drops for length (< 1000 bp).
  3. Assigns each to a curated reference label — accepted/tentative STs plus
     the provisional novel_N clusters — by 8-mer k-mer CONTAINMENT against
     every reference sequence, with top-5 consensus voting.

Why containment rather than Jaccard: a 600 bp fragment of an 1800 bp
reference shares at most ~1/3 of the union of their k-mer sets, so Jaccard
penalises short queries even on a perfect match. Containment
(1 - |Q ∩ R| / |Q|) measures how much of the QUERY is explained by the
reference, which is the right question for a fragment. Both strands are
tested and the better one kept.

Confidence tiers (containment distance to nearest reference):
  high ≤ 0.05   ·   medium ≤ 0.15   ·   low ≤ 0.30   ·   unmapped > 0.30
A "low" call additionally requires the top-5 nearest references to agree on
the label; disagreement is reported as ambiguous.

Usage (after a main-pipeline run has produced the reference FASTAs):
  python map_short_reads.py --email you@example.com \\
      --refs ./rRNA_pipeline_20260813/fasta \\
      --output-dir ./short_read_mapping

  # Re-analyse a cached download without hitting NCBI again:
  python map_short_reads.py --email you@example.com \\
      --refs ./rRNA_pipeline_20260813/fasta \\
      --output-dir ./short_read_mapping --skip-download
"""

import argparse
import importlib.util
import logging
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from Bio import SeqIO

# ── Import the main pipeline module for its query + metadata machinery ───────
_HERE = Path(__file__).resolve().parent
_PIPE_PATH = _HERE / "blastocystis_rRNA_pipeline_v5.py"
spec = importlib.util.spec_from_file_location("blasto_pipe", _PIPE_PATH)
pipe = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pipe)

K = 8
NBITS = 4 ** K  # 65,536 possible 8-mers
_BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
_COMP = str.maketrans("ACGT", "TGCA")


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--email", required=True, help="NCBI Entrez email")
    p.add_argument("--api-key", default="", help="NCBI API key")
    p.add_argument("--refs", required=True,
                   help="fasta/ dir of a main-pipeline run (*_filtered.fasta)")
    p.add_argument("--output-dir", default="./short_read_mapping")
    p.add_argument("--max-records", type=int, default=25000,
                   help="cap on the broad query")
    p.add_argument("--min-len", type=int, default=100,
                   help="shortest record worth k-mer profiling")
    p.add_argument("--short-max-len", type=int, default=999,
                   help="upper bound of the 'short' class (main pipeline "
                        "takes over from 1000 bp)")
    p.add_argument("--skip-download", action="store_true",
                   help="reuse <output-dir>/short_records.fasta + "
                        "short_metadata.csv from a previous run")
    p.add_argument("--thresh-high", type=float, default=0.05)
    p.add_argument("--thresh-medium", type=float, default=0.15)
    p.add_argument("--thresh-low", type=float, default=0.30)
    return p.parse_args()


def kmer_bitrow(seq: str) -> np.ndarray:
    """Indices of the distinct ACGT 8-mers present in seq."""
    s = seq.upper()
    idxs = set()
    for i in range(len(s) - K + 1):
        code = 0
        ok = True
        for ch in s[i:i + K]:
            b = _BASE.get(ch)
            if b is None:
                ok = False
                break
            code = (code << 2) | b
        if ok:
            idxs.add(code)
    return np.fromiter(idxs, dtype=np.int64, count=len(idxs))


def build_sparse(rows_of_indices) -> sparse.csr_matrix:
    indptr = [0]
    indices = []
    for r in rows_of_indices:
        indices.append(r)
        indptr.append(indptr[-1] + len(r))
    if indices:
        data = np.ones(indptr[-1], dtype=np.float32)
        indices = np.concatenate(indices) if indptr[-1] else np.array([], dtype=np.int64)
    else:
        data = np.array([], dtype=np.float32)
        indices = np.array([], dtype=np.int64)
    return sparse.csr_matrix((data, indices, np.array(indptr)),
                             shape=(len(indptr) - 1, NBITS))


def load_references(refs_dir: Path, log):
    """Load *_filtered.fasta; label = filename stem before '_filtered'.
    ST_unassigned is excluded — those are the singletons even the full-length
    analysis could not place, so they are not a meaningful mapping target."""
    ref_seqs, ref_labels = [], []
    for fa in sorted(refs_dir.glob("*_filtered.fasta")):
        label = fa.stem.replace("_filtered", "")
        if label == "ST_unassigned":
            continue
        for rec in SeqIO.parse(str(fa), "fasta"):
            ref_seqs.append(str(rec.seq))
            ref_labels.append(label)
    n_st = len({l for l in ref_labels if not l.startswith("novel_")})
    n_novel = len({l for l in ref_labels if l.startswith("novel_")})
    log.info(f"  References: {len(ref_seqs)} sequences, "
             f"{n_st} ST labels + {n_novel} novel clusters")
    return ref_seqs, ref_labels


def main():
    args = parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s  %(levelname)-8s  %(message)s",
        handlers=[logging.StreamHandler(),
                  logging.FileHandler(out / "mapping_log.txt")])
    log = logging.getLogger("short_map")

    pipe.Entrez.email = args.email
    if args.api_key:
        pipe.Entrez.api_key = args.api_key

    short_fa = out / "short_records.fasta"
    short_meta_csv = out / "short_metadata.csv"

    # ── 1. Acquire the short records ─────────────────────────────────────────
    if args.skip_download and short_fa.exists() and short_meta_csv.exists():
        log.info("[1] --skip-download: reusing cached short records")
        meta = pd.read_csv(short_meta_csv)
        seqs = {r.id: str(r.seq) for r in SeqIO.parse(str(short_fa), "fasta")}
        meta = meta[meta.accession.isin(seqs)]
    else:
        log.info("[1] Broad NCBI query (same query as the main pipeline)")
        records = pipe.fetch_records_by_query(
            pipe.broad_blastocystis_query(), args.max_records, log)
        n_total = len(records)
        shorts, too_short, full_len, too_long = [], 0, 0, 0
        for rec in records:
            L = len(rec.seq)
            if L < args.min_len:
                too_short += 1
            elif L <= args.short_max_len:
                shorts.append(rec)
            elif L <= 2500:
                full_len += 1
            else:
                too_long += 1
        log.info(f"  {n_total} records: {len(shorts)} short "
                 f"[{args.min_len}-{args.short_max_len} bp] to map, "
                 f"{full_len} full-length (main pipeline's territory), "
                 f"{too_short} < {args.min_len} bp, {too_long} > 2500 bp")

        rows = []
        for rec in shorts:
            md = pipe.extract_metadata_block(rec)
            md["accession"] = rec.id
            rows.append(md)
        meta = pd.DataFrame(rows)
        meta.to_csv(short_meta_csv, index=False)
        with open(short_fa, "w") as fh:
            for rec in shorts:
                fh.write(f">{rec.id}\n{rec.seq}\n")
        seqs = {rec.id: str(rec.seq) for rec in shorts}
        log.info(f"  Cached {len(seqs)} short records to {short_fa.name} "
                 f"(reuse with --skip-download)")

    # ── 2. Reference index ───────────────────────────────────────────────────
    log.info("[2] Building reference k-mer index")
    ref_seqs, ref_labels = load_references(Path(args.refs), log)
    R = build_sparse([kmer_bitrow(s) for s in ref_seqs])

    # ── 3. Containment assignment, both strands ──────────────────────────────
    log.info(f"[3] Assigning {len(meta)} short records by 8-mer containment")
    accs = meta.accession.tolist()
    fwd_rows, rc_rows, q_sizes = [], [], []
    for a in accs:
        s = seqs[a]
        f = kmer_bitrow(s)
        r = kmer_bitrow(s.translate(_COMP)[::-1])
        fwd_rows.append(f)
        rc_rows.append(r)
        q_sizes.append(max(len(f), 1))
    Qf = build_sparse(fwd_rows)
    Qr = build_sparse(rc_rows)

    inter = np.maximum((Qf @ R.T).toarray(), (Qr @ R.T).toarray())
    dist = 1.0 - inter / np.array(q_sizes, dtype=np.float32)[:, None]

    order5 = np.argsort(dist, axis=1)[:, :5]
    results = []
    for i, acc in enumerate(accs):
        top5_idx = order5[i]
        top5 = [(ref_labels[j], float(dist[i, j])) for j in top5_idx]
        best_label, d = top5[0]
        vote = Counter(l for l, _ in top5).most_common(1)[0][0]
        if d <= args.thresh_high:
            conf = "high"
        elif d <= args.thresh_medium:
            conf = "medium"
        elif d <= args.thresh_low:
            conf = "low" if vote == best_label else "ambiguous"
        else:
            conf, best_label = "unmapped", "unmapped"
        results.append({
            "accession": acc, "assigned": best_label,
            "containment_dist": round(d, 4), "confidence": conf,
            "top5": "; ".join(f"{l}:{dd:.3f}" for l, dd in top5),
        })

    res = pd.DataFrame(results)
    manifest = meta.merge(res, on="accession")
    manifest.to_csv(out / "manifest_short_reads.csv", index=False)

    # ── 4. Summaries ─────────────────────────────────────────────────────────
    log.info("[4] Summaries")
    mapped = manifest[~manifest.confidence.isin(["unmapped", "ambiguous"])]
    summ = (mapped.groupby("assigned")
            .agg(n=("accession", "size"),
                 median_dist=("containment_dist", "median"),
                 high=("confidence", lambda c: (c == "high").sum()),
                 top_hosts=("host", lambda h: "; ".join(
                     f"{k}:{v}" for k, v in Counter(h).most_common(3))),
                 top_countries=("country", lambda c: "; ".join(
                     f"{k}:{v}" for k, v in Counter(c).most_common(3))))
            .sort_values("n", ascending=False))
    summ.to_csv(out / "summary_by_label.csv")

    novel_hits = mapped[mapped.assigned.str.startswith("novel_")]
    novel_hits.to_csv(out / "novel_cluster_hits.csv", index=False)

    conf_counts = manifest.confidence.value_counts()
    log.info(f"  Mapped {len(mapped)}/{len(manifest)} "
             f"({len(mapped) / len(manifest) * 100:.1f}%) — " +
             "; ".join(f"{k}:{v}" for k, v in conf_counts.items()))
    log.info(f"  Barcode records hitting provisional novel clusters: "
             f"{len(novel_hits)} across "
             f"{novel_hits.assigned.nunique()} clusters")
    log.info("  Top 10 assigned labels:\n" +
             summ.head(10).to_string())

    # ── 5. Comparison chart: barcode world vs curated full-length world ──────
    log.info("[5] Comparison chart")
    plot_comparison(mapped, Path(args.refs), out, log)
    log.info(f"Done. Outputs in {out}/")


def plot_comparison(mapped, refs_dir, out, log):
    """Short-read ST distribution vs the curated full-length set's."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    full_counts = Counter()
    for fa in refs_dir.glob("*_filtered.fasta"):
        label = fa.stem.replace("_filtered", "")
        if label == "ST_unassigned":
            continue
        full_counts[label] = sum(1 for _ in SeqIO.parse(str(fa), "fasta"))

    short_counts = mapped.assigned.value_counts()
    top = short_counts.head(15).index.tolist()
    short_pct = [short_counts.get(l, 0) / short_counts.sum() * 100 for l in top]
    full_total = sum(full_counts.values())
    full_pct = [full_counts.get(l, 0) / full_total * 100 for l in top]

    fig, ax = plt.subplots(figsize=(8.5, 7), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    y = np.arange(len(top))
    h = 0.38
    ax.barh(y - h / 2 - 0.02, short_pct, height=h, color="#2a78d6",
            label="Barcode records (this mapping)", zorder=3)
    ax.barh(y + h / 2 + 0.02, full_pct, height=h, color="#eb6834",
            label="Curated full-length set", zorder=3)
    ax.set_yticks(y, top)
    ax.invert_yaxis()
    ax.set_xlabel("Share of records (%)", color="#52514e")
    ax.set_title("Subtype distribution: barcode surveillance vs curated "
                 "full-length references", color="#0b0b0b", loc="left",
                 fontsize=11)
    ax.legend(frameon=False, loc="lower right")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(colors="#52514e")
    ax.xaxis.grid(True, color="#e5e4e0", zorder=0)
    fig.tight_layout()
    fig.savefig(out / "st_distribution_short_vs_full.png", dpi=200)
    plt.close(fig)
    log.info(f"  ✔ {out / 'st_distribution_short_vs_full.png'}")


if __name__ == "__main__":
    main()
