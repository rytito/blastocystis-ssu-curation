# Blastocystis SSU rRNA curation pipeline

A single open-source Python pipeline that fetches, validates, places, clusters,
and reports on **every _Blastocystis_ SSU rRNA record on NCBI GenBank** — with
full provenance for every call it makes.

Presented at the 5th International _Blastocystis_ Conference (La Laguna,
Tenerife, September 2026): *"What the public database tells us — when we let
it: a curated re-reading of every Blastocystis SSU rRNA record on GenBank"*
(Raul Y. Tito Tadeo).

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

## Headline results (run of 2026-08-13)

The provenance manifest and per-ST summary of the latest run are in
[`results/run_20260813/`](results/run_20260813/).

- 18,074 records pulled by the broad query → **1,052** near-full-length,
  deduplicated, audited sequences across 62 ST labels (92.7 % known).
- Tight intra-ST clustering: median intra-distance 0.021 vs inter 0.121
  (Cliff's δ = −0.981) — STs are real, separable clusters.
- 18 provisional `novel_N` clusters; the tortoise cluster (7 records, six
  tortoise species, two countries: KT438705–KT438710 + EF209018) is
  consistent with a genuine chelonian lineage and is offered as a hypothesis
  for formal follow-up under the Stensvold & Clark criteria.

Numbers shift as new deposits land; the pipeline is designed to be re-run.
Between the May and August 2026 runs, 10 new records entered the curated set,
including the first direct ST42 record (Bos taurus, Portugal) and a new
2-record novel cluster from Eurasian beaver (*Castor fiber*, Czechia).

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

Run of 2026-08-13 ([`results/short_read_mapping_20260813/`](results/short_read_mapping_20260813/)):
16,784 barcode records, **92.8 % subtype-called** (12,010 high-confidence);
ST3/ST1/ST2 dominate as expected. Notably, **121 barcode records from
independent studies map into 11 of the 18 provisional novel clusters** —
including 7 high-confidence tortoise-lineage hits (four exact-match) from
accession series disjoint from the cluster's full-length references, and an
independent *Castor fiber* (Slovakia) hit on the beaver cluster first seen
in the August full-length run.

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

## Contributing

Curation cases, novel-cluster validations, and candidate accessions for the
reference table are very welcome — please open an issue.

## License

MIT — see [LICENSE](LICENSE).
