# Blastocystis SSU rRNA curation pipeline

A single open-source Python pipeline that fetches, validates, places, clusters,
and reports on **every _Blastocystis_ SSU rRNA record on NCBI GenBank** — with
full provenance for every call it makes.

## Why

_Blastocystis_ subtyping depends almost entirely on SSU rRNA. When a paper
reports "ST3 dominant in this cohort", the implicit chain of trust runs:
PCR → sequence → GenBank's nearest hit → that submitter's subtype call →
their reading of someone else's earlier call. That chain had never been
audited at scale. This pipeline applies a single transparent curation logic
to the whole database and reports what survives, what gets dropped, and why.

## What it does

1. **Acquire** — one broad NCBI Entrez query
   (`Blastocystis[Organism]` + 18S/SSU/rRNA title terms), no per-ST narrow
   queries. Paginated efetch with retry/backoff and resume support.
2. **Filter** — drop RefSeq-predicted transcripts, wrong `mol_type`, and
   barcode-region amplicons outside 1000–2500 bp (the ~600 bp Scicluna
   barcode lacks the variable regions needed for novel-ST detection).
3. **Assign ST** — two stages per record:
   (i) parse `/isolate`, `/clone`, `/strain`, `/note` for an explicit ST tag
   (high/medium/low confidence); (ii) for unlabelled records, 8-mer Jaccard
   placement against curated + high-confidence seeds with top-5 consensus.
4. **Deduplicate** — cross-ST MD5 deduplication catches identical sequences
   claimed under multiple STs.
5. **Cluster & report** — unplaceable records are k-mer-clustered into
   provisional `novel_N` groups; the pipeline emits a distance matrix
   (JC-corrected), hierarchical clustering, UMAP/PCoA ordination,
   per-ST summaries, host/geographic/temporal/technology distributions,
   and a full per-record manifest.

**Provenance:** every row of `manifest.csv` carries the NCBI query that found
it, its ST-confidence tier and validation reason, which GenBank field its host
came from, the source of its country, its MD5, and its retrieval timestamp —
everything needed to reproduce or contest the call.

## Install

```bash
conda env create -f environment.yml
conda activate blasto
```

External tools (installed by `environment.yml`): `mafft` (aligner),
`iqtree` (optional, for a publication-grade tree — suggested, not invoked).

## Run

```bash
# Full-scale acquisition — the configuration used for the conference talk
python blastocystis_rRNA_pipeline_v5.py \
    --email you@example.com \
    --max-broad 20000 \
    --jc-correct \
    --threads 4 \
    --seed 42

# Resume an interrupted run
python blastocystis_rRNA_pipeline_v5.py --email you@example.com --resume

# Re-analyse existing FASTAs without downloading
python blastocystis_rRNA_pipeline_v5.py --email you@example.com \
    --skip-download --jc-correct
```

An NCBI Entrez email is required; an `--api-key` raises the rate limit from
3 to 10 requests/s. The exact configuration of the talk run is preserved in
`config_talk_run_20260520.json` (every run writes its own `run_config.json`).

## Outputs

```
rRNA_pipeline/
├── fasta/manifest.csv        # per-record provenance manifest (the key output)
├── fasta/*.fasta             # per-ST and combined FASTAs
├── aligned/                  # MAFFT alignment
├── distances/                # pairwise distance matrix (JC-corrected)
├── clusters/cluster_assignments.csv
├── trees/                    # NJ tree (gated; use IQ-TREE for publication)
├── plots/                    # host, geography/time, technology, ordination,
│                             # novel clusters, interactive HTML report
├── summary_by_st.csv         # per-ST statistics
├── run_config.json           # exact arguments of this run
└── pipeline_log.txt
```

## Headline results (run of 2026-09-02)

The provenance manifest, per-ST summary, cluster assignments and full log of
the latest run are in [`results/run_20260902/`](results/run_20260902/); the
previous run is kept alongside in
[`results/run_20260813/`](results/run_20260813/) so the two can be diffed.

- 18,097 records pulled by the broad query -> **1,053** near-full-length,
  deduplicated, audited sequences across 61 labels (43 accepted or tentative
  subtypes covering 979 records).
- Tight intra-ST clustering: median intra-distance 0.021 vs inter 0.123
  (Cliff's delta = -0.976, over 478,731 pairs) - STs are real, separable clusters.
- 81 records match no accepted subtype: **17** provisional `novel_N` clusters
  (48 records) plus 26 `ST_unassigned` singletons. The tortoise cluster
  (7 records, six tortoise species, two countries: KT438705-KT438710 +
  EF209018) is consistent with a genuine chelonian lineage and is offered as a
  hypothesis for formal follow-up under the Stensvold & Clark criteria.

**`novel_N` labels are run-specific.** They are assigned in size order and
renumber between runs: the chelonian cluster is `novel_3` in the August run and
`novel_2` in September; the *Castor fiber* cluster `novel_10` -> `novel_9`.
Identify clusters by host/accession composition, never by number.

### What changed since 2026-08-13

Two things changed at once - a fresh GenBank screening and one methods change.
`compare_runs.py` and `compare_barcode.py` regenerate the full diff; the
narrative version is in [`CHANGES_20260902.md`](CHANGES_20260902.md).

- **One new record entered the set**, and it is the first African
  near-full-length sequence: `PZ873232.1`, ST2, Malawi (Zomba), 2024, ONT.
  Africa was a flat zero across 28 countries in August; it is now 1 record of
  910 with a known country, across 29.
- **A novel cluster dissolved.** With the enlarged ST2 seed set, two
  *Gorilla gorilla* records (`JX159034`, `JX159024`) and one Argentinian record
  (`MZ783086`) became placeable in ST2, so unplaced records fell 90 -> 81 and
  clusters 18 -> 17. Novel-cluster calls are sensitive to reference
  completeness; this is the same effect as the ST41->ST17 case.
- ST2 is the only subtype whose count changed (83 -> 87). Hosts, geography,
  technology composition and the separation statistics are otherwise stable.

## Companion analysis: barcode-region mapping

`map_short_reads.py` handles the ~16,000 records the main pipeline excludes
for length — the ~600 bp Scicluna barcode region that dominates surveillance
sequencing. It re-issues the same broad query, keeps the 100–999 bp records,
and assigns each to the curated reference labels (STs **and** the provisional
novel clusters) by 8-mer k-mer containment (1 − |Q∩R|/|Q|, both strands,
top-5 consensus — containment rather than Jaccard because a barcode fragment
shares at most a third of the k-mer union of a full-length reference even on
a perfect match).

```bash
python map_short_reads.py --email you@example.com \
    --refs ./rRNA_pipeline/fasta --output-dir ./short_read_mapping
```

### Separation check (added 2026-09-02)

Confidence tiers measure how *close* a fragment is to its nearest reference,
not whether that reference is distinguishable from the runner-up. Because ST10
was split into ST42-44 on a ~2 % full-length threshold, and a ~600 bp window
cannot resolve that, a fragment can sit inside the "high" band of two subtypes
at once: in the August run **100 % of ST44 calls and 84 % of ST43 calls** had a
sibling in their own top-5, a quarter of them exact ties.

`--margin` (default 0.01) now requires every call to beat the nearest
*differently-labelled* reference by that margin. Failing calls collapse to an
agreed complex label (`ST10-complex = {ST10, ST42, ST43, ST44}`) or, where no
complex is defined, are marked `ambiguous`. A `separation_note` column records
which sibling caused each collapse and by what gap. Thanks to Eleni Gentekaki
for raising this.

Run of 2026-09-02 ([`results/short_read_mapping_20260902/`](results/short_read_mapping_20260902/)):
16,806 barcode records, **90.3 % subtype-called** (11,711 high-confidence),
down from 92.8 % in August because of the separation check. 178 records are now
reported as `ST10-complex`; ST44 drops 174 -> 57 and ST43 58 -> 48, while ST42
barely moves (403 -> 397), consistent with its 3.2-3.9 % separation from the
other three.

**68 barcode records from independent studies map into 6 of the 17 provisional
novel clusters** (was 121 across 11). The corroboration that fell away was
colliding with an accepted subtype - novel_5 vs ST4, novel_8 vs ST3, novel_4 vs
ST12 - so it was never evidence of novelty. The chelonian lineage is untouched:
12 records, 7 at high confidence, before and after.

## What it is not

- Not a replacement for expert curation — a confidently wrong submitter tag
  still passes stage-I validation.
- The curated set itself excludes ~16,000 Scicluna-region barcode records by
  design — they lack the variable regions needed to define novel STs. They
  are covered separately by `map_short_reads.py` (above), which calls them
  against the curated references instead.
- `novel_N` labels are clustering labels, not subtype proposals.
- The 8-mer Jaccard placement is calibrated, not optimal — plug in EPA-ng or
  pplacer for formal phylogenetic placement.


## License

MIT — see [LICENSE](LICENSE).
