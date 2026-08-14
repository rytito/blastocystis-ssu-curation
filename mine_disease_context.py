"""
Mine disease context from Blastocystis GenBank records
=======================================================

Question: do the ~18,000 Blastocystis SSU rRNA deposits carry enough
disease metadata to probe subtype-disease associations in humans —
colorectal cancer in particular?

For every record returned by the main pipeline's broad query, this script
extracts the free-text fields where deposit context lives:

  - DEFINITION line
  - source qualifiers: /note, /isolation_source, /host
  - titles of the REFERENCE blocks (the studies behind each deposit)

and classifies each record against a disease keyword taxonomy. Subtype
labels are joined from the two existing manifests (curated full-length set
+ barcode mapping), so every record that received an ST call anywhere in
the project can be cross-tabulated as ST x disease context.

IMPORTANT CAVEATS (also printed into the report):
  - Record-level qualifiers (/note, /isolation_source) describe the actual
    sample; REFERENCE titles describe the STUDY. A record matched only via
    its study title (context="study") may be a control sample from a
    disease-focused cohort. Both levels are kept, separately.
  - GenBank metadata are voluntary and sparse; absence of a disease term
    is not evidence of health. Records tagged healthy/control are the only
    usable comparison group.
  - This is deposit-metadata mining, not epidemiology. Output is
    hypothesis-generating only.

Usage:
  python mine_disease_context.py --email you@example.com \
      --curated ./rRNA_pipeline/fasta/manifest.csv \
      --barcode ./short_read_mapping/manifest_short_reads.csv \
      --output-dir ./disease_context
"""

import argparse
import importlib.util
import logging
import re
from collections import Counter
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "blasto_pipe", _HERE / "blastocystis_rRNA_pipeline_v5.py")
pipe = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pipe)

# ── Disease taxonomy: category -> regex (case-insensitive) ───────────────────
DISEASE_PATTERNS = {
    "colorectal_cancer": r"colorectal|colon\s+cancer|rectal\s+cancer|\bCRC\b|"
                         r"adenoma|colonic\s+neoplas|bowel\s+cancer|"
                         r"intestinal\s+polyp",
    "other_cancer":      r"\bcancer\b|carcinoma|leukemi|lymphoma|malignan|"
                         r"tumou?r|oncolog|chemotherap",
    "ibs":               r"irritable\s+bowel|\bIBS\b",
    "ibd":               r"inflammatory\s+bowel|\bIBD\b|crohn|ulcerative\s+colitis",
    "diarrhea_gi":       r"diarrh|gastroenteritis|dysenter|gastrointestinal\s+"
                         r"symptom|abdominal\s+pain|digestive\s+(disorder|symptom)|"
                         r"\bvomit|flatulen|bloating|dyspep",
    "urticaria_skin":    r"urticaria|\brash\b|pruritus|angioedema|dermatolog",
    "immunocompromised": r"\bHIV\b|\bAIDS\b|immunocompromis|immunosuppress|"
                         r"immunodefic|transplant|h(a)?emodialysis|dialysis",
    "anemia_nutrition":  r"an(a)?emi[ac]|malnutrition|malnourish|stunting|iron\s+"
                         r"(status|deficien)",
    "healthy_control":   r"healthy|asymptomatic|control\s+(group|subject|individual)|"
                         r"non[-\s]?symptomatic|without\s+(gastrointestinal\s+)?symptom",
    "symptomatic_unspec": r"symptomatic\s+(patient|subject|individual|carrier)|"
                          r"\bpatient(s)?\b",
}
CATEGORY_ORDER = list(DISEASE_PATTERNS)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--email", required=True)
    p.add_argument("--api-key", default="")
    p.add_argument("--curated", required=True,
                   help="manifest.csv of the curated full-length run")
    p.add_argument("--barcode", required=True,
                   help="manifest_short_reads.csv of the barcode mapping")
    p.add_argument("--output-dir", default="./disease_context")
    p.add_argument("--max-records", type=int, default=25000)
    p.add_argument("--skip-download", action="store_true",
                   help="reuse <output-dir>/record_context.csv")
    return p.parse_args()


def classify(text: str) -> dict:
    hits = {}
    for cat, pat in DISEASE_PATTERNS.items():
        m = re.search(pat, text, flags=re.IGNORECASE)
        if m:
            hits[cat] = m.group(0)
    # "patient(s)" alone is weak; drop symptomatic_unspec if a stronger
    # category matched
    if "symptomatic_unspec" in hits and len(hits) > 1:
        del hits["symptomatic_unspec"]
    # colorectal beats generic cancer
    if "colorectal_cancer" in hits and "other_cancer" in hits:
        del hits["other_cancer"]
    return hits


def main():
    args = parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s  %(levelname)-8s  %(message)s",
                        handlers=[logging.StreamHandler(),
                                  logging.FileHandler(out / "mining_log.txt")])
    log = logging.getLogger("disease")

    ctx_csv = out / "record_context.csv"
    if args.skip_download and ctx_csv.exists():
        log.info("[1] Reusing cached record context")
        ctx = pd.read_csv(ctx_csv)
    else:
        pipe.Entrez.email = args.email
        if args.api_key:
            pipe.Entrez.api_key = args.api_key
        log.info("[1] Fetching all records (broad query, full GenBank format)")
        records = pipe.fetch_records_by_query(
            pipe.broad_blastocystis_query(), args.max_records, log)
        rows = []
        for rec in records:
            note = isol = host = ""
            for feat in rec.features:
                if feat.type == "source":
                    note = " ".join(str(n) for n in feat.qualifiers.get("note", []))
                    isol = feat.qualifiers.get("isolation_source", [""])[0]
                    host = feat.qualifiers.get("host", [""])[0]
                    break
            titles = "; ".join(
                (r.title or "") for r in rec.annotations.get("references", []))
            host_norm, _ = pipe.coalesce_host_from_gb(rec)
            rows.append({
                "accession": rec.id,
                "length": len(rec.seq),
                "host": pipe.normalise_host(host_norm),
                "definition": rec.description,
                "note": note, "isolation_source": isol, "host_raw": host,
                "ref_titles": titles,
            })
        ctx = pd.DataFrame(rows)
        ctx.to_csv(ctx_csv, index=False)
        log.info(f"  cached {len(ctx)} records -> {ctx_csv.name}")

    # ── 2. Classify: record-level vs study-level evidence ────────────────────
    log.info("[2] Classifying disease context")
    rec_cats, rec_evid, study_cats, study_evid = [], [], [], []
    for _, r in ctx.iterrows():
        rec_text = " | ".join(str(x) for x in
                              [r.definition, r.note, r.isolation_source, r.host_raw]
                              if pd.notna(x))
        st_text = str(r.ref_titles) if pd.notna(r.ref_titles) else ""
        rh = classify(rec_text)
        sh = classify(st_text)
        rec_cats.append(";".join(rh) if rh else "")
        rec_evid.append(";".join(f"{k}={v}" for k, v in rh.items()))
        study_cats.append(";".join(sh) if sh else "")
        study_evid.append(";".join(f"{k}={v}" for k, v in sh.items()))
    ctx["record_categories"] = rec_cats
    ctx["record_evidence"] = rec_evid
    ctx["study_categories"] = study_cats
    ctx["study_evidence"] = study_evid

    # ── 3. Join subtype calls from both manifests ────────────────────────────
    log.info("[3] Joining subtype assignments")
    cur = pd.read_csv(args.curated)[["accession", "subtype"]]
    cur["st_source"] = "curated_fulllength"
    bar = pd.read_csv(args.barcode)[["accession", "assigned", "confidence"]]
    bar = bar[~bar.confidence.isin(["unmapped", "ambiguous"])]
    bar = bar.rename(columns={"assigned": "subtype"})[["accession", "subtype"]]
    bar["st_source"] = "barcode_mapping"
    st = pd.concat([cur, bar]).drop_duplicates("accession")
    ctx = ctx.merge(st, on="accession", how="left")
    ctx.to_csv(out / "record_context_classified.csv", index=False)

    # ── 4. Cross-tabulations (human records with an ST call) ────────────────
    log.info("[4] Cross-tabulating ST x disease (human records)")
    hum = ctx[(ctx.host == "Homo sapiens") & ctx.subtype.notna()].copy()
    log.info(f"  human records with ST call: {len(hum)}")

    def explode_counts(frame, col):
        f = frame[frame[col] != ""].copy()
        f = f.assign(cat=f[col].str.split(";")).explode("cat")
        return f.reset_index(drop=True)

    for level in ["record", "study"]:
        col = f"{level}_categories"
        e = explode_counts(hum, col)
        n_annot = hum[hum[col] != ""].shape[0]
        log.info(f"  {level}-level: {n_annot}/{len(hum)} human ST-called "
                 f"records carry any disease annotation "
                 f"({n_annot / max(len(hum), 1) * 100:.1f}%)")
        xt = pd.crosstab(e.subtype, e.cat)
        xt.to_csv(out / f"st_by_disease_{level}_level.csv")
        log.info(f"  categories ({level}): " +
                 "; ".join(f"{k}:{v}" for k, v in
                           Counter(e.cat).most_common()))

    # CRC drill-down, both levels
    for level in ["record", "study"]:
        col = f"{level}_categories"
        crc = hum[hum[col].str.contains("colorectal_cancer", na=False)]
        crc.to_csv(out / f"crc_records_{level}_level.csv", index=False)
        log.info(f"  CRC-context ({level}-level): {len(crc)} records | STs: " +
                 ("; ".join(f"{k}:{v}" for k, v in
                            Counter(crc.subtype).most_common()) or "none"))

    log.info(f"Done. Outputs in {out}/")


if __name__ == "__main__":
    main()
