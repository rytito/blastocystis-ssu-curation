"""
Blastocystis SSU rRNA Pipeline v5 — Broad-first acquisition at full scale
==========================================================================

Comprehensive acquisition of all available Blastocystis SSU rRNA sequences
from NCBI, with content-based ST assignment, hierarchical host coalescing,
and full provenance tracking.

What changed in v5 vs v4:
  • BROAD-FIRST acquisition. v4 ran 40 per-ST narrow queries (capped at
    200 each) plus a single capped broad query, missing most of the
    ~3,000 records currently in NCBI for Blastocystis 18S rRNA. v5 makes
    the broad genus-wide query the primary source and partitions records
    by ST via content validation + placement. Default cap raised to
    --max-broad 10000 (effectively uncapped for current corpus size).
  • Robust paginated efetch. Large pulls (>500 records) use 100-per-batch
    fetches with explicit retry, exponential backoff, and progress
    checkpointing so a network hiccup does not lose the run.
  • Two-stage ST assignment for every record:
      (i) parse /isolate, /clone, /strain, /note for explicit ST tag
          → high/medium/low confidence (same as v4)
      (ii) for records with no ST tag → 8-mer Jaccard placement against
          curated reference seeds + previously high-confidence records
      Records that get no usable assignment in either stage are kept
      as "ST_unassigned" with their metadata, in a separate FASTA
      that downstream phylogenetic analysis can place explicitly.
  • Optional per-ST narrow query as a coverage CHECK rather than the
    primary acquisition path (--narrow-coverage-check, default off).
    Catches any records the broad query missed (rare but possible).
  • --resume support reads any existing manifest.csv and avoids
    re-fetching accessions already processed.

Everything else from v4 carries forward unchanged:
  • Hierarchical host coalescing (/host → /isolation_source → /note → organism)
  • Patched sequencing-technology detection (IonTorrent, BGI/MGI, etc.)
  • Cross-ST MD5 deduplication
  • Vectorised p-distance, JC correction
  • Distance-from-medoid novel-variant flagging
  • Country / geo_loc_name dual extraction
  • Natural ST sort order in all plots

External tools (in PATH):
  mafft     — recommended aligner
  muscle    — fallback (v3 and v5 CLIs supported)
  iqtree    — for publication-grade tree (suggested, not invoked)

Usage:
  # Full default acquisition (~3,000 records, ~10-20 min)
  python blastocystis_rRNA_pipeline_v5.py --email you@example.com --jc-correct

  # Resume an interrupted run
  python blastocystis_rRNA_pipeline_v5.py --email you@example.com --resume

  # Add per-ST narrow query as a coverage cross-check
  python blastocystis_rRNA_pipeline_v5.py --email you@example.com \\
      --narrow-coverage-check

  # Skip acquisition entirely (re-analyse existing FASTAs)
  python blastocystis_rRNA_pipeline_v5.py --email you@example.com \\
      --skip-download --jc-correct
"""

import argparse
import csv
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path

# ── Dependency checks ─────────────────────────────────────────────────────────
MISSING = []
try:
    import numpy as np
except ImportError:
    MISSING.append("numpy")
try:
    import pandas as pd
except ImportError:
    MISSING.append("pandas")
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
except ImportError:
    MISSING.append("matplotlib seaborn")
try:
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import linkage, fcluster
except ImportError:
    MISSING.append("scipy")
try:
    from Bio import Entrez, SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    MISSING.append("biopython")
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError:
    MISSING.append("plotly")
try:
    from tqdm import tqdm
except ImportError:
    MISSING.append("tqdm")

if MISSING:
    sys.exit(
        f"[ERROR] Missing packages: {', '.join(MISSING)}\n"
        f"Run: pip install {' '.join(MISSING)}\n"
        f"  or: conda install -c conda-forge {' '.join(MISSING)}"
    )

try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  REFERENCE DATA                                                          ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# ── Blastocystis subtype registry ─────────────────────────────────────────────
SUBTYPE_REGISTRY = {
    "ST1":  {"status": "accepted",   "human": True,  "note": "Most prevalent globally; IBS association"},
    "ST2":  {"status": "accepted",   "human": True,  "note": "Common in humans; urticaria association"},
    "ST3":  {"status": "accepted",   "human": True,  "note": "Most prevalent globally; IBS-D association"},
    "ST4":  {"status": "accepted",   "human": True,  "note": "Europe-dominant; possible immunoprotective role"},
    "ST5":  {"status": "accepted",   "human": True,  "note": "Predominantly pigs; zoonotic in humans"},
    "ST6":  {"status": "accepted",   "human": True,  "note": "Primarily birds/poultry"},
    "ST7":  {"status": "accepted",   "human": True,  "note": "Primarily birds/poultry"},
    "ST8":  {"status": "accepted",   "human": True,  "note": "Non-human primates; rare in humans"},
    "ST9":  {"status": "accepted",   "human": True,  "note": "Human-specific; rare"},
    "ST10": {"status": "accepted",   "human": True,  "note": "Artiodactyls; rare human reports; now split ST42-44"},
    "ST11": {"status": "accepted",   "human": False, "note": "Non-human primates"},
    "ST12": {"status": "accepted",   "human": True,  "note": "Non-human primates; occasional in humans"},
    "ST13": {"status": "accepted",   "human": True,  "note": "Rare; occasional human reports"},
    "ST14": {"status": "accepted",   "human": True,  "note": "Cattle/artiodactyls; rare in humans"},
    "ST15": {"status": "accepted",   "human": False, "note": "Non-human primates"},
    "ST16": {"status": "accepted",   "human": True,  "note": "Rare; occasional human reports"},
    "ST17": {"status": "accepted",   "human": False, "note": "Rodents; cattle"},
    "ST18": {"status": "contested",  "human": False, "note": "Proposed invalid — likely chimeric"},
    "ST19": {"status": "contested",  "human": False, "note": "Proposed invalid — likely chimeric"},
    "ST20": {"status": "contested",  "human": False, "note": "Proposed invalid — likely chimeric"},
    "ST21": {"status": "accepted",   "human": False, "note": "Avian; accepted after full-length validation"},
    "ST22": {"status": "contested",  "human": False, "note": "Proposed invalid — likely chimeric"},
    "ST23": {"status": "accepted",   "human": True,  "note": "Reported in humans (Colombia, Thailand)"},
    "ST24": {"status": "tentative",  "human": False, "note": "Cattle (USA)"},
    "ST25": {"status": "tentative",  "human": False, "note": "Cattle"},
    "ST26": {"status": "accepted",   "human": False, "note": "Cattle (USA); dominant in several US states"},
    "ST27": {"status": "accepted",   "human": False, "note": "Full-length MinION validated"},
    "ST28": {"status": "accepted",   "human": False, "note": "Full-length MinION validated"},
    "ST29": {"status": "accepted",   "human": False, "note": "Full-length MinION validated"},
    "ST30": {"status": "tentative",  "human": False, "note": "Limited data"},
    "ST31": {"status": "tentative",  "human": False, "note": "Limited data"},
    "ST32": {"status": "tentative",  "human": False, "note": "Limited data"},
    "ST33": {"status": "tentative",  "human": False, "note": "Limited data"},
    "ST34": {"status": "tentative",  "human": False, "note": "Limited data"},
    "ST35": {"status": "accepted",   "human": True,  "note": "Human isolate 168; Nanopore full-length validated"},
    "ST36": {"status": "accepted",   "human": False, "note": "Wild bat (S. lilium)"},
    "ST37": {"status": "accepted",   "human": False, "note": "Wild rodent (Heteromyidae)"},
    "ST38": {"status": "accepted",   "human": False, "note": "European water vole"},
    "ST39": {"status": "tentative",  "human": False, "note": "Recently proposed"},
    "ST40": {"status": "tentative",  "human": False, "note": "Recently proposed; Algeria study"},
    "ST41": {"status": "tentative",  "human": True,  "note": "Reported in humans (China 2024)"},
    "ST42": {"status": "accepted",   "human": False, "note": "ST10 split; Santin et al. 2024"},
    "ST43": {"status": "accepted",   "human": False, "note": "ST10 split; Santin et al. 2024"},
    "ST44": {"status": "accepted",   "human": False, "note": "ST10 split; Santin et al. 2024"},
    "NMAST1": {"status": "avian_nmast", "human": False, "note": "Non-mammalian avian ST"},
    "NMAST2": {"status": "avian_nmast", "human": False, "note": "Non-mammalian avian ST"},
    "NMAST3": {"status": "avian_nmast", "human": False, "note": "Non-mammalian avian ST"},
    "NMAST4": {"status": "avian_nmast", "human": False, "note": "Non-mammalian avian ST"},
    "NMAST5": {"status": "avian_nmast", "human": False, "note": "Non-mammalian avian ST"},
    "NMAST6": {"status": "avian_nmast", "human": False, "note": "Non-mammalian avian ST"},
    "NMAST7": {"status": "avian_nmast", "human": False, "note": "Non-mammalian avian ST"},
    "NMAST8": {"status": "avian_nmast", "human": False, "note": "Non-mammalian avian ST"},
}

SUBTYPES_DEFAULT = [st for st in SUBTYPE_REGISTRY if not st.startswith("NMAST")]


# ── Literature-curated reference accessions (pass-2a seeds) ──────────────────
# Extend conservatively — every entry should trace to a peer-reviewed
# ST assignment. Used as labelled seeds for phylogenetic placement.
REFERENCE_ST_MAP = {
    "AB070987": "ST1",   # Arisue 2003, HJ96-1 strain (B. hominis from human)
    "U51151":   "ST1",   # Silberman 1996, NAND II canonical ST1 reference
    # Add more as you confirm them from publications. Suggested next:
    # "AF408425": "ST3",   # Stensvold 2007 ST3 reference
    # "AY244620": "ST4",   # Yoshikawa 2004 ST4 reference
    # "DQ366344": "ST6",   # Avian ST6 reference
}


# ── Subtype colour palette ────────────────────────────────────────────────────
import colorsys

def _make_palette(n):
    colors = []
    for i in range(n):
        h = i / n
        s = 0.75 if i % 2 == 0 else 0.55
        v = 0.85 if (i // 2) % 2 == 0 else 0.65
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        colors.append(f"#{int(r*255):02X}{int(g*255):02X}{int(b*255):02X}")
    return colors

_ALL_ST_KEYS = [k for k in SUBTYPE_REGISTRY if not k.startswith("NMAST")]
SUBTYPE_COLORS = dict(zip(_ALL_ST_KEYS, _make_palette(len(_ALL_ST_KEYS))))
SUBTYPE_COLORS.update({
    "ST1": "#E41A1C", "ST2": "#377EB8", "ST3": "#4DAF4A",
    "ST4": "#984EA3", "ST5": "#FF7F00", "ST6": "#A65628",
    "ST7": "#F781BF", "ST8": "#999999", "ST9": "#66C2A5",
})
# Provisional novel-subtype clusters get a bright, attention-getting palette.
# Distinct from accepted STs so they're obvious in figures.
NOVEL_COLORS = ["#FF1493", "#00CED1", "#FFD700", "#9370DB",
                "#00FA9A", "#FF6347", "#1E90FF", "#FF8C00"]
for i, c in enumerate(NOVEL_COLORS, start=1):
    SUBTYPE_COLORS[f"novel_{i}"] = c
SUBTYPE_COLORS["ST_unassigned"] = "#444444"  # dark grey

STATUS_MARKERS = {
    "accepted":    {"marker": "o", "alpha": 0.9, "edgecolor": "none"},
    "tentative":   {"marker": "^", "alpha": 0.7, "edgecolor": "gray"},
    "contested":   {"marker": "x", "alpha": 0.5, "edgecolor": "red"},
    "avian_nmast": {"marker": "D", "alpha": 0.6, "edgecolor": "none"},
    "novel":       {"marker": "*", "alpha": 0.95, "edgecolor": "black"},
    "unassigned":  {"marker": ".", "alpha": 0.5,  "edgecolor": "none"},
}


# ── Sequencing technology detection ──────────────────────────────────────────
TECH_PATTERNS = [
    ("Sanger", [
        r"sanger", r"capillar[y]?\s+sequenc",
        r"abi\s*(?:prism|3[0-9]{2,3}(?:xl)?)",
        r"\b3500xl?\b", r"\b3730xl?\b",
        r"\b377\b\s*(?:dna\s*sequenc)?",
        r"dye.?terminator", r"\bbigdye\b",
        r"conventional\s+sequenc", r"direct\s+sequenc",
        r"chain.?termination",
    ]),
    ("Illumina", [
        r"illumina",
        r"\bmiseq\b", r"\bhiseq\b", r"\bnextseq\b", r"\bnovaseq\b",
        r"\biseq\b", r"\bminiseq\b",
        r"\bgenome\s+analyz?er\b", r"\bgaii+\b",
        r"short.?read\s+sequenc", r"paired.?end\s+sequenc",
    ]),
    ("ONT", [
        r"oxford\s+nanopore", r"\bont\b", r"\bnanopore\b",
        r"\bminion\b", r"\bgridion\b", r"\bpromethion\b", r"\bflongle\b",
        r"long.?read.*nanopore",
    ]),
    ("PacBio", [
        r"pacbio", r"pacific\s+biosciences", r"\bsmrt\b",
        r"\bsequel\s*ii?e?\b", r"\brs\s*ii\b", r"\brsii\b", r"\brevio\b",
    ]),
    ("Ion Torrent", [
        r"\bion\s*torrent\b", r"\biontorrent\b",
        r"\bion\s*pgm\b", r"\bion\s*proton\b",
        r"\bion\s*(?:gene)?studio\b", r"\bion\s*s5\b", r"\bion\s*chef\b",
        r"thermo.*ion", r"life\s*tech.*ion",
    ]),
    ("BGI/MGI", [
        r"\bbgi\s*seq", r"\bbgiseq\b",
        r"\bmgi\s*seq", r"\bmgiseq\b",
        r"\bdnb\s*seq", r"\bdnbseq\b",
        r"\bbgi[\s-]*tech", r"\bcomplete\s+genomics\b",
    ]),
    ("Element AVITI", [
        r"\belement\s+aviti\b", r"\baviti\b", r"\belement\s+bioscience",
    ]),
    ("Singular G4", [
        r"\bsingular\s+g4\b", r"\bsingular\s+genomics\b",
    ]),
    ("Ultima", [r"\bultima\s+genomics\b", r"\bug100\b"]),
    ("454/Pyrosequencing", [
        r"\b454\b", r"pyrosequenc", r"roche\s+454",
        r"gs.?flx", r"gs.?junior",
        r"genome\s+sequencer\s+(?:flx|junior)",
    ]),
    ("SOLiD", [
        r"\bsolid\s*(?:[345]\d{2}|system|sequenc)",
        r"applied\s+biosystems\s+solid",
    ]),
    ("Helicos", [
        r"\bhelicos\b", r"\bheliscope\b", r"true\s+single\s+molecule",
    ]),
]

STRUCTURED_COMMENT_KEY_FRAGMENTS = [
    "sequencing technology", "sequencing-technology", "sequencing_technology",
    "sequencing platform", "sequencing-platform",
    "assembly method", "instrument", "platform",
    "sequencing instrument", "sequencing kit",
    "library construction", "library preparation",
]

TECH_COLORS = {
    "Sanger": "#377EB8", "Illumina": "#4DAF4A", "ONT": "#FF7F00",
    "PacBio": "#984EA3", "Ion Torrent": "#A65628",
    "BGI/MGI": "#17BECF", "Element AVITI": "#BCBD22",
    "Singular G4": "#9467BD", "Ultima": "#8C564B",
    "454/Pyrosequencing": "#F781BF", "SOLiD": "#999999",
    "Helicos": "#7F7F7F", "Unknown": "#DDDDDD",
}


# ── Host coalescing data ──────────────────────────────────────────────────────
ORGANISM_HOST_MAP = {
    "blastocystis hominis":  "Homo sapiens",
    "blastocystis ratti":    "Rattus norvegicus",
    "blastocystis galli":    "Gallus gallus",
    "blastocystis lapemi":   "Lapemis hardwickii",
    "blastocystis pythoni":  "Python reticulatus",
    "blastocystis anatis":   "Anas platyrhynchos",
    "blastocystis anseri":   "Anser anser",
}

GENERIC_HOSTS = {
    "", "unknown", "n/a", "na", "not applicable", "none", "missing",
    "not collected", "not available", "blastocystis sp.",
    "blastocystis hominis",
}

HOST_NORMALISATION = {
    # Humans
    "human": "Homo sapiens", "humans": "Homo sapiens",
    "homo sapien": "Homo sapiens", "homo sapiens": "Homo sapiens",
    # Cattle
    "cattle": "Bos taurus", "cow": "Bos taurus", "bovine": "Bos taurus",
    "bos taurus": "Bos taurus", "calf": "Bos taurus", "beef cattle": "Bos taurus",
    "dairy cattle": "Bos taurus",
    # Pigs
    "pig": "Sus scrofa", "swine": "Sus scrofa",
    "domestic pig": "Sus scrofa domesticus", "sus scrofa": "Sus scrofa",
    "sus scrofa domesticus": "Sus scrofa domesticus", "piglet": "Sus scrofa domesticus",
    # Sheep
    "sheep": "Ovis aries", "domestic sheep": "Ovis aries",
    "ovis aries": "Ovis aries", "lamb": "Ovis aries", "ewe": "Ovis aries",
    # Goats
    "goat": "Capra hircus", "capra hircus": "Capra hircus",
    # Chicken / poultry
    "chicken": "Gallus gallus", "gallus gallus": "Gallus gallus",
    "fowl": "Gallus gallus", "broiler": "Gallus gallus",
    "domestic chicken": "Gallus gallus",
    # Camels
    "camel": "Camelus dromedarius", "dromedary": "Camelus dromedarius",
    "camelus dromedarius": "Camelus dromedarius",
    # Equids
    "horse": "Equus caballus", "donkey": "Equus asinus",
    # Pets
    "dog": "Canis lupus familiaris", "cat": "Felis catus",
    # Primates
    "gibbon": "Hylobates lar",
    "monkey": "Primates sp.", "macaque": "Macaca sp.",
    "chimpanzee": "Pan troglodytes", "gorilla": "Gorilla gorilla",
    "orangutan": "Pongo pygmaeus",
    # Rodents
    "rat": "Rattus sp.", "mouse": "Mus musculus",
    "rattus norvegicus": "Rattus norvegicus",
    # Other commonly seen
    "muskox": "Ovibos moschatus", "ovibos moschatus": "Ovibos moschatus",
    "ctenodactylus gundi": "Ctenodactylus gundi",
    "yak": "Bos grunniens", "buffalo": "Bubalus bubalis",
    "water buffalo": "Bubalus bubalis",
    "alpaca": "Vicugna pacos", "llama": "Lama glama",
    "reindeer": "Rangifer tarandus", "elk": "Cervus canadensis",
    "deer": "Cervidae sp.",
    "duck": "Anas platyrhynchos", "goose": "Anser anser",
    "turkey": "Meleagris gallopavo",
    "bat": "Chiroptera sp.",
}

NOTE_HOST_PATTERNS = [
    re.compile(r"\bisolate(?:d)?\s+from\s+(?:the\s+)?(?:a\s+)?([a-zA-Z][a-zA-Z\s]+?)\b",
               re.IGNORECASE),
    re.compile(r"\bfrom\s+(?:a\s+)?human\b", re.IGNORECASE),
    re.compile(r"\bhuman\s+(?:fecal|stool|isolate)\b", re.IGNORECASE),
]


# ── ST tag patterns (for validation step) ─────────────────────────────────────
ST_TAG_PATTERNS = [
    re.compile(r"\bsubtype\s*(\d+)\b", re.IGNORECASE),
    re.compile(r"\bST\s*(\d+)\b"),
    re.compile(r"\b(?:Blastocystis\s+sp\.?\s*)?ST(\d+)[-_]", re.IGNORECASE),
]


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  METADATA EXTRACTION FROM GENBANK                                        ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def infer_sequencing_technology(gb_record) -> str:
    """Extract sequencing technology from a GenBank SeqRecord."""
    search_texts = []

    sc = gb_record.annotations.get("structured_comment", {})
    if isinstance(sc, dict):
        for block_name, block_dict in sc.items():
            if not isinstance(block_dict, dict):
                continue
            for key, val in block_dict.items():
                key_lower = key.lower()
                if any(frag in key_lower for frag in STRUCTURED_COMMENT_KEY_FRAGMENTS):
                    search_texts.append(str(val).lower())

    comment = gb_record.annotations.get("comment", "")
    if comment:
        search_texts.append(str(comment).lower())

    for feat in gb_record.features:
        if feat.type in ("source", "misc_feature", "gene"):
            for key in ("note", "db_xref", "mol_type", "isolation_source",
                        "PCR_primers"):
                for val in feat.qualifiers.get(key, []):
                    search_texts.append(str(val).lower())

    if gb_record.description:
        search_texts.append(gb_record.description.lower())
    if gb_record.name:
        search_texts.append(gb_record.name.lower())

    combined = " ".join(search_texts)
    for tech_name, patterns in TECH_PATTERNS:
        for pat in patterns:
            if re.search(pat, combined, re.IGNORECASE):
                return tech_name
    return "Unknown"


def normalise_host(s: str) -> str:
    """Light host-string normalisation. Returns the input unchanged if no rule applies."""
    if not isinstance(s, str):
        return "Unknown"
    s_clean = s.strip().lower()
    if s_clean in GENERIC_HOSTS:
        return "Unknown"
    return HOST_NORMALISATION.get(s_clean, s.strip())


def coalesce_host_from_gb(gb_rec) -> tuple[str, str]:
    """
    Return (host_string, source_of_host) where source_of_host is one of:
      "host_qualifier", "isolation_source", "note", "organism", "unknown"

    Priority order:
      1. /host qualifier, if not generic
      2. /isolation_source if it names an animal (not a sample type)
      3. /note text matching "isolated from X" patterns
      4. Organism binomial (Blastocystis hominis -> human, etc.)
    """
    host_q = isol_q = note_q = ""
    organism = (gb_rec.annotations.get("organism") or "").strip().lower()

    for feat in gb_rec.features:
        if feat.type != "source":
            continue
        host_q = (feat.qualifiers.get("host",             [""])[0] or "").strip()
        isol_q = (feat.qualifiers.get("isolation_source", [""])[0] or "").strip()
        notes  = feat.qualifiers.get("note", [])
        note_q = " ".join(str(n) for n in notes).strip()
        break

    # 1. /host qualifier
    if host_q and host_q.strip().lower() not in GENERIC_HOSTS:
        return host_q, "host_qualifier"

    # 2. /isolation_source — accept if it names an animal
    if isol_q and isol_q.strip().lower() not in GENERIC_HOSTS:
        s_low = isol_q.strip().lower()
        looks_like_animal = (
            s_low in HOST_NORMALISATION
            or re.search(r"\b(human|cattle|cow|pig|swine|sheep|goat|chicken|"
                         r"camel|horse|dog|cat|bird|monkey|gibbon|rat|mouse|"
                         r"muskox|gundi|reindeer|bovine|porcine|caprine|ovine|"
                         r"poultry|fowl|deer|buffalo|yak|alpaca|llama|"
                         r"chimpanzee|gorilla|orangutan|macaque|bat|hominid)\b",
                         s_low)
            or re.match(r"^[a-z]+\s+[a-z]+$", s_low)
        )
        # Exclude sample types (the most common reason /isolation_source is
        # not the host)
        looks_like_sample = re.search(
            r"\b(feces|faeces|stool|soil|water|sewage|sludge|manure|gut|"
            r"intestine|colon|cecum|environment|culture|laboratory|"
            r"fecal\s+sample|fecal\s+matter)\b", s_low
        )
        if looks_like_animal and not looks_like_sample:
            return isol_q, "isolation_source"

    # 3. /note text
    if note_q:
        for pat in NOTE_HOST_PATTERNS:
            m = pat.search(note_q)
            if m:
                grp = m.group(1) if m.groups() else "Human"
                grp = grp.strip()
                if grp.lower() not in GENERIC_HOSTS:
                    return grp.title(), "note"

    # 4. Organism binomial
    for org_key, host in ORGANISM_HOST_MAP.items():
        if org_key in organism:
            return host, "organism"

    return "Unknown", "unknown"


def extract_metadata_block(gb_rec) -> dict:
    """Extract a complete metadata block from a single GenBank record.

    Note on country: INSDC renamed the /country qualifier to /geo_loc_name
    in June 2024 (mandatory from Dec 2024). Older records use /country;
    newer ones use /geo_loc_name. We check both and prefer whichever is
    populated. See https://ncbiinsights.ncbi.nlm.nih.gov/2023/12/14/update-genbank-qualifier/
    """
    country = date_ = isol = "Unknown"
    for feat in gb_rec.features:
        if feat.type == "source":
            # Try /geo_loc_name first (current standard), then /country (legacy)
            geo = feat.qualifiers.get("geo_loc_name", [""])[0].strip()
            old = feat.qualifiers.get("country",      [""])[0].strip()
            country = geo or old or "Unknown"
            date_   = feat.qualifiers.get("collection_date", ["Unknown"])[0]
            isol    = feat.qualifiers.get("isolation_source",["Unknown"])[0]
            break

    host_raw, host_source = coalesce_host_from_gb(gb_rec)
    host = normalise_host(host_raw)
    tech = infer_sequencing_technology(gb_rec)
    md5 = hashlib.md5(str(gb_rec.seq).upper().encode()).hexdigest()

    return {
        "host":             host,
        "host_source":      host_source,
        "country":          country,
        "collection_date":  date_,
        "isolation_source": isol,
        "seq_tech":         tech,
        "md5":              md5,
        "length":           len(gb_rec.seq),
    }


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  ST CLAIM VALIDATION                                                     ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def extract_declared_sts(gb_rec) -> set:
    """
    Return the set of STs explicitly declared in a GenBank record's source
    feature qualifiers (isolate, clone, strain, note) plus the first ~120
    characters of the description. Excludes the abstract and publication
    title — those routinely co-mention other STs.
    """
    declared = set()
    bits = []

    for feat in gb_rec.features:
        if feat.type != "source":
            continue
        for key in ("isolate", "clone", "strain", "note", "specific_host"):
            for val in feat.qualifiers.get(key, []):
                bits.append(str(val))

    bits.append(gb_rec.name or "")
    bits.append((gb_rec.description or "")[:120])

    blob = " ".join(bits)
    for pat in ST_TAG_PATTERNS:
        for m in pat.findall(blob):
            try:
                declared.add(f"ST{int(m)}")
            except ValueError:
                pass
    return declared


def validate_st_claim(gb_rec, queried_st: str) -> tuple[bool, str, str]:
    """
    Decide whether to keep a fetched record for the queried ST.

    Returns (keep, confidence, reason).
      keep=True   → use this record for queried_st
      confidence ∈ {"high", "medium", "low", "rejected"}
        high   : source declares exactly queried_st
        medium : source declares queried_st alongside others (mixed)
        low    : no ST declared anywhere; fall back to query trust
        rejected: source declares only STs ≠ queried_st
    """
    if queried_st.startswith("NMAST"):
        return True, "low", "NMAST — query trusted"

    declared = extract_declared_sts(gb_rec)
    if not declared:
        return True, "low", "no ST declared in source feature; query trusted"

    if queried_st in declared:
        if len(declared) == 1:
            return True, "high", f"source declares {queried_st} only"
        others = sorted(declared - {queried_st})
        return True, "medium", f"source declares {queried_st} alongside {others}"

    return False, "rejected", f"source declares {sorted(declared)}, not {queried_st}"


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  ARGUMENT PARSING & LOGGING                                              ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def parse_args():
    p = argparse.ArgumentParser(
        description="Blastocystis SSU rRNA Pipeline v5 — broad-first comprehensive acquisition",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # NCBI / I/O
    g = p.add_argument_group("Input / output")
    g.add_argument("--email",      type=str, default="",   help="NCBI Entrez email (required for download)")
    g.add_argument("--api-key",    type=str, default="",   help="NCBI API key (10 req/s vs 3 req/s)")
    g.add_argument("--output-dir", type=str, default="./rRNA_pipeline")
    g.add_argument("--skip-download", action="store_true",
                   help="Reanalyse existing FASTAs in <output-dir>/fasta/")
    g.add_argument("--resume", action="store_true",
                   help="Skip re-fetching accessions already in manifest.csv")

    # Subtype selection (only used if --narrow-coverage-check is on)
    g = p.add_argument_group("Subtype selection (narrow coverage check only)")
    g.add_argument("--subtypes", type=str, default="",
                   help="Comma-separated STs (default: all from registry)")
    g.add_argument("--include-nmast",     action="store_true")
    g.add_argument("--include-contested", action="store_true")
    g.add_argument("--status-filter",     type=str, default="",
                   help="Comma-separated status values to keep")

    # Length filter
    g = p.add_argument_group("Length filter")
    g.add_argument("--min-len", type=int, default=1000)
    g.add_argument("--max-len", type=int, default=2500)

    # Acquisition strategy — v5 broad-first design
    g = p.add_argument_group("Acquisition strategy (v5 broad-first)")
    g.add_argument("--max-broad", type=int, default=10000,
                   help="Cap on broad query records (default: 10000, effectively uncapped)")
    g.add_argument("--fetch-batch-size", type=int, default=100,
                   help="Records per efetch call (larger = faster, more memory)")
    g.add_argument("--fetch-retries", type=int, default=3,
                   help="Retries per failed fetch batch with exponential backoff")
    g.add_argument("--narrow-coverage-check", action="store_true",
                   help="Also run per-ST narrow queries as a coverage cross-check")
    g.add_argument("--max-per-st", type=int, default=200,
                   help="Cap for narrow per-ST query (only used with --narrow-coverage-check)")
    g.add_argument("--max-placement-dist", type=float, default=0.30,
                   help="k-mer Jaccard cutoff for ST placement of unlabelled records")
    g.add_argument("--no-st-validation", action="store_true",
                   help="Skip ST-claim validation (trust NCBI metadata blindly)")
    g.add_argument("--no-cross-st-dedup", action="store_true",
                   help="Don't deduplicate identical sequences across STs")
    g.add_argument("--keep-unassigned", action="store_true", default=True,
                   help="Keep records that fail ST assignment in a separate ST_unassigned.fasta")

    # Novel-subtype clustering of unassigned records
    g.add_argument("--cluster-novel-subtypes", action="store_true", default=True,
                   help="Cluster ST_unassigned records into provisional novel_N groups")
    g.add_argument("--novel-cluster-cutoff", type=float, default=0.20,
                   help="k-mer Jaccard cutoff for novel-subtype clustering (tighter than placement)")
    g.add_argument("--novel-min-size", type=int, default=2,
                   help="Minimum cluster size to receive a novel_N label (singletons stay ST_unassigned)")

    # Aligner
    g = p.add_argument_group("Alignment")
    g.add_argument("--aligner", type=str, default="auto",
                   choices=["mafft", "muscle", "kmer", "auto"])
    g.add_argument("--threads", type=int, default=4)

    # Distances & clustering
    g = p.add_argument_group("Distances & clustering")
    g.add_argument("--jc-correct",     action="store_true",
                   help="Apply Jukes-Cantor correction to p-distances")
    g.add_argument("--thresh-cluster", type=float, default=0.05)
    g.add_argument("--thresh-subtype", type=float, default=0.025)
    g.add_argument("--thresh-variant", type=float, default=0.01)
    g.add_argument("--novel-distance-multiplier", type=float, default=3.0,
                   help="Flag novel variants at >N× median intra-distance from medoid")

    # Misc
    g = p.add_argument_group("Miscellaneous")
    g.add_argument("--seed", type=int, default=42)
    g.add_argument("--verbose", action="store_true")

    return p.parse_args()


def setup_logging(out_dir: str, verbose: bool = False):
    log_path = os.path.join(out_dir, "pipeline_log.txt")
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    return logging.getLogger("pipeline")


def init_entrez(email: str, api_key: str, log):
    if not email:
        email = input("NCBI Entrez email (required): ").strip()
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
        log.info("NCBI API key set — 10 req/sec")
    else:
        log.info("No API key — 3 req/sec (use --api-key to speed up)")


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 1A — PASS 1: NARROW PER-ST QUERY                                   ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def build_narrow_query(subtype: str) -> str:
    """Per-ST title-scoped query. No broad [All Fields] subtype clause."""
    num = re.search(r"\d+", subtype)
    n = num.group() if num else subtype

    if subtype.startswith("NMAST"):
        nmast_n = subtype.replace("NMAST", "")
        clause = (
            f'"NMAST{nmast_n}"[Title] OR '
            f'"non-mammalian avian ST{nmast_n}"[Title] OR '
            f'"avian subtype {nmast_n}"[Title]'
        )
    else:
        clause = (
            f'"subtype {n}"[Title] OR "ST{n}"[Title] OR "ST {n}"[Title]'
        )

    return (
        f'Blastocystis[Organism] AND ({clause}) AND ('
        f'"18S"[Title] OR "SSU"[Title] OR "small subunit"[Title] OR '
        f'"ribosomal RNA"[Title] OR "rRNA"[Title] OR "rDNA"[Title] OR '
        f'"18S ribosomal"[All Fields] OR "SSU rRNA"[All Fields]'
        f') NOT ("whole genome"[Title] OR "complete genome"[Title] OR '
        f'"WGS"[Title] OR "shotgun"[Title])'
    )


def _entrez_sleep_interval() -> float:
    """NCBI rate-limit-friendly interval. 10 req/s with API key, 3 req/s without."""
    return 0.11 if getattr(Entrez, "api_key", None) else 0.4


def _retry_fetch(fetch_fn, retries: int, log, label: str = ""):
    """
    Run fetch_fn() with up to `retries` attempts on failure, with
    exponential backoff. fetch_fn must be a no-arg callable returning the
    handle's parsed result (caller decides what that is).
    Returns the result, or None if all attempts failed.
    """
    last_exc = None
    for attempt in range(retries + 1):
        try:
            return fetch_fn()
        except Exception as e:
            last_exc = e
            if attempt < retries:
                wait = 2 ** attempt   # 1s, 2s, 4s, 8s
                log.warning(f"  {label}: attempt {attempt + 1} failed ({e}); "
                            f"retrying in {wait}s...")
                time.sleep(wait)
    log.error(f"  {label}: all {retries + 1} attempts failed: {last_exc}")
    return None


def fetch_records_by_query(query: str, max_seqs: int, log,
                            batch_size: int = 100, retries: int = 3) -> list:
    """
    Run esearch + paginated efetch and return SeqRecords as GenBank objects.
    v5 differences vs v4:
      - default batch_size=100 (was 20) — ~5× faster on large pulls
      - per-batch retry with exponential backoff
      - reports total found vs total fetched at end
      - progress bar via tqdm for runs >200 records
    """
    sleep = _entrez_sleep_interval()

    # esearch with retry
    def _do_search():
        h = Entrez.esearch(db="nucleotide", term=query,
                           retmax=max_seqs, usehistory="y")
        result = Entrez.read(h); h.close()
        return result

    result = _retry_fetch(_do_search, retries, log, "esearch")
    if result is None:
        return []
    ids = result.get("IdList", [])
    webenv = result.get("WebEnv")
    query_key = result.get("QueryKey")
    total = int(result.get("Count", 0))

    if total > max_seqs:
        log.info(f"  {total} matches total — capped at --max-broad {max_seqs}; "
                 f"raise this if you want full coverage.")
    elif total == max_seqs:
        log.warning(f"  Fetched {total} which equals --max-broad cap; "
                    f"there may be more records — re-run with a higher cap.")
    if not ids:
        return []

    log.info(f"  Fetching {len(ids)} records in batches of {batch_size}...")

    records = []
    n_batches = (len(ids) + batch_size - 1) // batch_size
    use_pbar = len(ids) > 200

    batch_iter = range(0, len(ids), batch_size)
    if use_pbar:
        batch_iter = tqdm(batch_iter, total=n_batches, desc="  efetch",
                          unit="batch", leave=False)

    for start in batch_iter:
        def _do_batch():
            h = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                              retstart=start, retmax=batch_size,
                              webenv=webenv, query_key=query_key)
            batch = list(SeqIO.parse(h, "genbank"))
            h.close()
            return batch

        batch = _retry_fetch(_do_batch, retries, log,
                             f"efetch batch [{start}:{start + batch_size}]")
        if batch:
            records.extend(batch)
        time.sleep(sleep)

    log.info(f"  Retrieved {len(records)}/{len(ids)} records "
             f"({total} total matches in NCBI).")
    return records


def fetch_accessions(accessions: list, log, batch_size: int = 100,
                     retries: int = 3) -> dict:
    """
    Fetch GenBank records by accession ID list. Returns dict keyed by record.id.
    v5 uses batch_size=100 by default and adds per-batch retry.
    """
    sleep = _entrez_sleep_interval()
    accs = list(dict.fromkeys(accessions))
    records: dict = {}

    if not accs:
        return records

    log.info(f"  Fetching {len(accs)} accessions in batches of {batch_size}...")
    n_batches = (len(accs) + batch_size - 1) // batch_size
    use_pbar = len(accs) > 200

    batch_iter = range(0, len(accs), batch_size)
    if use_pbar:
        batch_iter = tqdm(batch_iter, total=n_batches, desc="  efetch",
                          unit="batch", leave=False)

    for start in batch_iter:
        batch_accs = accs[start:start + batch_size]

        def _do_batch():
            h = Entrez.efetch(db="nucleotide", id=",".join(batch_accs),
                              rettype="gb", retmode="text")
            batch = list(SeqIO.parse(h, "genbank"))
            h.close()
            return batch

        batch = _retry_fetch(_do_batch, retries, log,
                             f"accession batch [{start}:{start + batch_size}]")
        if batch:
            for rec in batch:
                records[rec.id] = rec
        time.sleep(sleep)

    log.info(f"  Retrieved {len(records)}/{len(accs)} records.")
    return records


def pass1_fetch_per_st(subtype: str, max_seqs: int, validate_st: bool, log):
    """Pass 1: per-ST narrow query with optional content validation."""
    query = build_narrow_query(subtype)
    records = fetch_records_by_query(query, max_seqs, log)
    if not records:
        log.info(f"  {subtype}: 0 sequences found")
        return [], [], query

    log.info(f"  {subtype}: {len(records)} sequences fetched — validating...")

    kept_records = []
    manifest_rows = []
    decision_counts = Counter()
    retrieved_at = datetime.now(timezone.utc).isoformat()

    for gb_rec in records:
        if validate_st:
            keep, confidence, reason = validate_st_claim(gb_rec, subtype)
            decision_counts[confidence] += 1
            if not keep:
                continue
        else:
            confidence, reason = "unchecked", "validation disabled"

        meta = extract_metadata_block(gb_rec)
        fasta_rec = SeqRecord(
            gb_rec.seq, id=gb_rec.id,
            description=(
                f"{gb_rec.description} [subtype={subtype}] "
                f"[host={meta['host']}] [country={meta['country']}] "
                f"[date={meta['collection_date']}] "
                f"[isolation_source={meta['isolation_source']}] "
                f"[seq_tech={meta['seq_tech']}] "
                f"[confidence={confidence}] [host_source={meta['host_source']}]"
            )
        )
        kept_records.append(fasta_rec)
        manifest_rows.append({
            "accession":         gb_rec.id,
            "subtype":           subtype,
            "length":            meta["length"],
            "md5":               meta["md5"],
            "host":              meta["host"],
            "host_source":       meta["host_source"],
            "country":           meta["country"],
            "collection_date":   meta["collection_date"],
            "isolation_source":  meta["isolation_source"],
            "seq_tech":          meta["seq_tech"],
            "st_confidence":     confidence,
            "st_validation":     reason,
            "acquisition_pass":  "pass1_narrow",
            "retrieved_at_utc":  retrieved_at,
            "ncbi_query":        query,
        })

    if validate_st and decision_counts:
        rejected = decision_counts.get("rejected", 0)
        summary = ", ".join(f"{k}:{v}" for k, v in decision_counts.items())
        if rejected:
            log.info(f"    ST-validation: {rejected} rejected as co-mention "
                     f"contamination ({summary})")
        else:
            log.info(f"    ST-validation: all records pass ({summary})")

    log.info(f"  {subtype}: {len(kept_records)} kept after validation")
    return kept_records, manifest_rows, query


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 1B — PASS 2: BROAD QUERY + PHYLOGENETIC PLACEMENT                  ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def broad_blastocystis_query() -> str:
    """
    Genus-wide Blastocystis SSU rRNA query, no ST filtering.

    v5 design note: this query intentionally allows 18S/SSU/rRNA tokens
    anywhere in the record ([All Fields]) rather than restricting to
    [Title]. NCBI's default expansion of the user-facing query
    'Blastocystis AND 18S rRNA' returns ~3,000 records this way as of
    May 2026; v4's title-restricted version returned far fewer.
    The trade-off — a few non-SSU-rRNA records may slip through — is
    handled downstream by the length filter (1000-2500 bp typical for
    near-full-length SSU rRNA) and by sequence-based ST validation.
    """
    return (
        'Blastocystis[Organism] AND ('
        '"18S"[All Fields] OR "SSU"[All Fields] OR '
        '"small subunit"[All Fields] OR '
        '"18S ribosomal"[All Fields] OR "SSU rRNA"[All Fields] OR '
        '"ribosomal RNA"[All Fields] OR "rRNA"[All Fields] OR '
        '"rDNA"[All Fields]'
        ') NOT ('
        '"whole genome"[Title] OR "complete genome"[Title] OR '
        '"WGS"[Title] OR "shotgun"[Title] OR '
        '"mitochondrial"[Title] OR "mitochondrion"[Title]'
        ')'
    )


def kmer_set(seq, k: int = 8) -> set:
    s = str(seq).upper().replace("-", "").replace("N", "")
    return set(s[i:i + k] for i in range(len(s) - k + 1))


def jaccard_dist(s1: set, s2: set) -> float:
    if not s1 or not s2:
        return 1.0
    return 1.0 - len(s1 & s2) / len(s1 | s2)


def place_against_seeds(unknown_records, seed_records, seed_st_map, log,
                        max_dist: float = 0.30) -> dict:
    """
    Assign each unknown record to its nearest-seed ST by 8-mer Jaccard
    distance. Top-5 nearest seeds vote for consensus.

    Returns {accession: (assigned_st, distance, confidence, reason)}.
    """
    log.info(f"  Indexing {len(seed_records)} seeds...")
    seed_kmers = {rec.id: kmer_set(rec.seq) for rec in seed_records}
    seed_ids = list(seed_kmers.keys())

    log.info(f"  Placing {len(unknown_records)} records against seeds...")
    placements = {}
    for rec in tqdm(unknown_records, desc="  Placement", leave=False):
        q_kmers = kmer_set(rec.seq)
        if not q_kmers:
            placements[rec.id] = ("Unplaced", 1.0, "rejected", "empty k-mer set")
            continue

        # Distance to every seed
        dists = [(jaccard_dist(q_kmers, seed_kmers[sid]), sid) for sid in seed_ids]
        dists.sort()
        best_dist, best_seed = dists[0]

        # Top-5 consensus
        top5_sts = [seed_st_map.get(sid, "Unknown") for _, sid in dists[:5]]
        top5_top = Counter(top5_sts).most_common(1)[0]
        assigned_st = seed_st_map.get(best_seed, "Unknown")

        # Determine confidence
        if top5_top[0] != assigned_st:
            confidence = "low"
            reason = (f"placement_low (nearest={assigned_st} d={best_dist:.3f}; "
                      f"top5={top5_sts})")
        elif best_dist > max_dist:
            confidence = "rejected"
            reason = f"placement_rejected (d={best_dist:.3f} > {max_dist})"
        elif best_dist > 0.15:
            confidence = "low"
            reason = f"placement_low (d={best_dist:.3f})"
        elif best_dist > 0.08:
            confidence = "medium"
            reason = f"placement_medium (d={best_dist:.3f})"
        else:
            confidence = "high"
            reason = f"placement_high (d={best_dist:.3f})"

        if confidence == "rejected":
            placements[rec.id] = ("Unplaced", best_dist, confidence, reason)
        else:
            placements[rec.id] = (assigned_st, best_dist, confidence, reason)

    return placements


# ── Novel-subtype clustering for unplaced records ────────────────────────────
def cluster_novel_subtypes(unplaced_records, log,
                            k: int = 8,
                            cluster_distance_cutoff: float = 0.20,
                            min_cluster_size: int = 2,
                            min_cluster_share: float = 0.6) -> dict:
    """
    Group records that the placement step could not assign to any known ST
    into provisional "novel_N" groups based on their own k-mer similarity.

    Method: pairwise 8-mer Jaccard distance, average-linkage hierarchical
    clustering, cut the dendrogram at `cluster_distance_cutoff` (default
    0.20 — tighter than the placement cutoff of 0.30, since we want
    coherent novel groups, not promiscuous lumping).

    Clusters smaller than `min_cluster_size` (default 2 — i.e. singletons)
    are NOT assigned a novel ST and remain `ST_unassigned`. This avoids
    giving every poorly-sequenced singleton its own "novel" ST.

    Numbering: clusters are sorted by size (largest first) and named
    novel_1, novel_2, ...; this makes novel_1 the most populous candidate
    group, which is the one most worth following up with a formal
    phylogenetic analysis.

    Returns: {accession: novel_label} for records that received a novel
    assignment. Records left as `ST_unassigned` are not in the returned
    dict.

    Important caveat (must be in the report / abstract): "novel_N" is a
    provisional cluster label, NOT a proposal of a new ST. Formal ST
    description requires the criteria laid out by Stensvold & Clark
    (~99% identity from a full-length SSU rRNA reference) plus PCR
    primers validated across multiple isolates. The novel_N labels are
    a hypothesis-generating tool, not a taxonomic claim.
    """
    if len(unplaced_records) < min_cluster_size:
        log.info(f"  Only {len(unplaced_records)} unplaced records — "
                 f"too few for novel clustering.")
        return {}

    log.info(f"  Novel clustering: {len(unplaced_records)} unplaced records "
             f"at k-mer Jaccard cutoff {cluster_distance_cutoff}")

    # Compute pairwise distances
    kmers = [kmer_set(rec.seq, k=k) for rec in unplaced_records]
    n = len(unplaced_records)
    # condensed distance matrix for scipy linkage
    condensed = []
    for i in range(n):
        for j in range(i + 1, n):
            condensed.append(jaccard_dist(kmers[i], kmers[j]))
    condensed = np.asarray(condensed)

    # Edge case: all records identical (distance 0 throughout)
    if condensed.size == 0:
        return {}

    # Average linkage on the condensed matrix
    Z = linkage(condensed, method="average")
    cluster_ids = fcluster(Z, t=cluster_distance_cutoff, criterion="distance")

    # Group records by cluster
    by_cluster = defaultdict(list)
    for i, cid in enumerate(cluster_ids):
        by_cluster[int(cid)].append(unplaced_records[i])

    # Keep only clusters with at least min_cluster_size members
    qualifying = [(cid, members) for cid, members in by_cluster.items()
                  if len(members) >= min_cluster_size]

    # Sort by size descending (largest cluster = novel_1)
    qualifying.sort(key=lambda x: -len(x[1]))

    # Compute median intra-cluster distance for each — a useful sanity stat
    novel_label_to_records: dict = {}
    novel_assignments: dict = {}
    for rank, (cid, members) in enumerate(qualifying, start=1):
        label = f"novel_{rank}"
        novel_label_to_records[label] = members
        intra = []
        member_ids = {r.id for r in members}
        for i, ri in enumerate(unplaced_records):
            if ri.id not in member_ids:
                continue
            for j, rj in enumerate(unplaced_records):
                if j <= i or rj.id not in member_ids:
                    continue
                intra.append(jaccard_dist(kmers[i], kmers[j]))
        med = float(np.median(intra)) if intra else 0.0
        log.info(f"    {label}: n={len(members)}  "
                 f"median intra={med:.3f}  "
                 f"example={members[0].id}")
        for r in members:
            novel_assignments[r.id] = label

    n_singletons = sum(1 for cid, members in by_cluster.items()
                       if len(members) < min_cluster_size)
    log.info(f"  Novel clustering result: "
             f"{len(qualifying)} provisional novel groups, "
             f"{n_singletons} singletons remain ST_unassigned.")

    return novel_assignments


def pass2_acquire(pass1_records, pass1_manifest, args, log) -> tuple:
    """Pass 2: literature-curated references + broad query + placement."""
    log.info("\n[STEP 1b] PASS 2 — broad acquisition + phylogenetic placement")

    pass1_accs = {m["accession"] for m in pass1_manifest}
    pass1_accs_bare = {a.split(".")[0] for a in pass1_accs}
    pass1_md5s = {m["md5"] for m in pass1_manifest}

    # ── 2a. Literature-curated references ───────────────────────────────────
    log.info("  Pass 2a — literature-curated reference accessions")
    ref_to_fetch = [a for a in REFERENCE_ST_MAP
                    if a not in pass1_accs_bare and a not in pass1_accs]
    log.info(f"    {len(REFERENCE_ST_MAP)} known refs, "
             f"{len(ref_to_fetch)} not in pass-1.")
    ref_records = fetch_accessions(ref_to_fetch, log) if ref_to_fetch else {}

    # ── 2b. Broad query ─────────────────────────────────────────────────────
    log.info("  Pass 2b — broad genus-wide query")
    broad_query = broad_blastocystis_query()
    try:
        h = Entrez.esearch(db="nucleotide", term=broad_query,
                           retmax=args.max_broad, usehistory="y")
        result = Entrez.read(h); h.close()
        broad_ids = result.get("IdList", [])
        log.info(f"    Broad query: {result.get('Count', 0)} matches; "
                 f"got {len(broad_ids)} IDs")
    except Exception as e:
        log.warning(f"    Broad query failed: {e}")
        broad_ids = []

    # Exclude what we already have
    broad_to_fetch = [
        i for i in broad_ids
        if i not in pass1_accs and i.split(".")[0] not in pass1_accs_bare
        and i.split(".")[0] not in REFERENCE_ST_MAP
        and i not in ref_records
    ]
    log.info(f"    {len(broad_to_fetch)} new accessions to fetch.")
    broad_records = fetch_accessions(broad_to_fetch, log) if broad_to_fetch else {}

    # ── Length-filter both batches ──────────────────────────────────────────
    def _len_filter(d, label):
        kept = {k: r for k, r in d.items()
                if args.min_len <= len(r.seq) <= args.max_len}
        log.info(f"    Length filter on {label}: {len(d)} → {len(kept)}")
        return kept

    ref_records   = _len_filter(ref_records,   "reference set")
    broad_records = _len_filter(broad_records, "broad query")

    # ── 2c. Build seed set ──────────────────────────────────────────────────
    log.info("  Pass 2c — building placement seeds from pass-1")
    seed_records = list(pass1_records)
    seed_st_map = {m["accession"]: m["subtype"] for m in pass1_manifest}
    # Add references to seed pool
    for acc, rec in ref_records.items():
        st = REFERENCE_ST_MAP.get(acc) or REFERENCE_ST_MAP.get(rec.id.split(".")[0])
        if st:
            seed_records.append(rec)
            seed_st_map[rec.id] = st
    log.info(f"    Seeds: {len(seed_records)} sequences spanning "
             f"{len(set(seed_st_map.values()))} STs.")

    # ── 2d. Place broad records against seeds ───────────────────────────────
    placements = {}
    if broad_records and seed_records:
        log.info("  Pass 2d — phylogenetic placement")
        placements = place_against_seeds(
            list(broad_records.values()), seed_records, seed_st_map, log,
            max_dist=args.max_placement_dist,
        )
        conf_counts = Counter()
        st_counts = Counter()
        for acc, (st, _, conf, _) in placements.items():
            conf_counts[conf] += 1
            if conf != "rejected":
                st_counts[st] += 1
        log.info("    Placement confidence breakdown:")
        for c, n in conf_counts.most_common():
            log.info(f"      {c:<10}: {n}")
        log.info("    Top placement STs:")
        for st, n in st_counts.most_common(8):
            log.info(f"      {st:<6}: {n}")

    # ── 2e. Build records + manifest rows for additions ─────────────────────
    log.info("  Pass 2e — merging pass-2 additions")
    retrieved_at = datetime.now(timezone.utc).isoformat()
    new_records = []
    new_manifest = []

    # Reference additions
    for acc, rec in ref_records.items():
        st = REFERENCE_ST_MAP.get(acc) or REFERENCE_ST_MAP.get(rec.id.split(".")[0])
        if not st:
            continue
        meta = extract_metadata_block(rec)
        if meta["md5"] in pass1_md5s:
            continue
        fasta_rec = SeqRecord(
            rec.seq, id=rec.id,
            description=(
                f"{rec.description} [subtype={st}] "
                f"[host={meta['host']}] [country={meta['country']}] "
                f"[date={meta['collection_date']}] "
                f"[isolation_source={meta['isolation_source']}] "
                f"[seq_tech={meta['seq_tech']}] "
                f"[confidence=curated_reference] [host_source={meta['host_source']}]"
            )
        )
        new_records.append(fasta_rec)
        new_manifest.append({
            "accession":         rec.id,
            "subtype":           st,
            "length":            meta["length"],
            "md5":               meta["md5"],
            "host":              meta["host"],
            "host_source":       meta["host_source"],
            "country":           meta["country"],
            "collection_date":   meta["collection_date"],
            "isolation_source":  meta["isolation_source"],
            "seq_tech":          meta["seq_tech"],
            "st_confidence":     "curated_reference",
            "st_validation":     f"literature-curated reference (mapped to {st})",
            "acquisition_pass":  "pass2_reference",
            "retrieved_at_utc":  retrieved_at,
            "ncbi_query":        "REFERENCE_LOOKUP",
        })
        pass1_md5s.add(meta["md5"])

    # Broad-query placements
    for acc, (st, dist, conf, reason) in placements.items():
        if st == "Unplaced":
            continue
        rec = broad_records.get(acc)
        if rec is None:
            continue
        meta = extract_metadata_block(rec)
        if meta["md5"] in pass1_md5s:
            continue
        fasta_rec = SeqRecord(
            rec.seq, id=rec.id,
            description=(
                f"{rec.description} [subtype={st}] "
                f"[host={meta['host']}] [country={meta['country']}] "
                f"[date={meta['collection_date']}] "
                f"[isolation_source={meta['isolation_source']}] "
                f"[seq_tech={meta['seq_tech']}] "
                f"[confidence={conf}] [host_source={meta['host_source']}]"
            )
        )
        new_records.append(fasta_rec)
        new_manifest.append({
            "accession":         rec.id,
            "subtype":           st,
            "length":            meta["length"],
            "md5":               meta["md5"],
            "host":              meta["host"],
            "host_source":       meta["host_source"],
            "country":           meta["country"],
            "collection_date":   meta["collection_date"],
            "isolation_source":  meta["isolation_source"],
            "seq_tech":          meta["seq_tech"],
            "st_confidence":     conf,
            "st_validation":     reason,
            "acquisition_pass":  "pass2_placement",
            "retrieved_at_utc":  retrieved_at,
            "ncbi_query":        broad_query,
        })
        pass1_md5s.add(meta["md5"])

    log.info(f"  Pass 2 added {len(new_records)} unique sequences "
             f"(after MD5 dedup vs pass-1).")
    return new_records, new_manifest


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 1C — CROSS-ST MD5 DEDUPLICATION                                    ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def deduplicate_within_st(records, log):
    """Per-ST exact-sequence dedup (gap- and case-insensitive)."""
    seen = {}
    for r in records:
        key = str(r.seq).upper().replace("-", "")
        if key not in seen:
            seen[key] = r
    if len(seen) < len(records):
        log.info(f"    Within-ST dedup: {len(records)} → {len(seen)}")
    return list(seen.values())


def deduplicate_across_sts(records, manifest_rows, log):
    """
    Global MD5 dedup across STs. Keeps the first occurrence of each unique
    sequence and reports which ST pairs had to share. Defense in depth
    against any cross-ST contamination that slipped past ST validation.
    """
    md5_by_acc = {m["accession"]: m["md5"] for m in manifest_rows}
    st_by_acc  = {m["accession"]: m["subtype"] for m in manifest_rows}

    seq_to_keep = {}
    seq_to_sts  = defaultdict(set)
    for r in records:
        md5 = md5_by_acc.get(r.id)
        if md5 is None:
            continue
        seq_to_sts[md5].add(st_by_acc.get(r.id, "?"))
        if md5 not in seq_to_keep:
            seq_to_keep[md5] = r

    duplicate_md5s = {md5: sts for md5, sts in seq_to_sts.items() if len(sts) > 1}
    if duplicate_md5s:
        log.warning(f"  {len(duplicate_md5s)} sequences shared across STs:")
        pair_counts = Counter()
        for sts in duplicate_md5s.values():
            sts = sorted(sts)
            for i in range(len(sts)):
                for j in range(i + 1, len(sts)):
                    pair_counts[(sts[i], sts[j])] += 1
        for (a, b), c in pair_counts.most_common(10):
            log.warning(f"    {a} ↔ {b}: {c} shared")

    deduped = list(seq_to_keep.values())
    log.info(f"  Cross-ST dedup: {len(records)} → {len(deduped)}")
    return deduped


def filter_by_length(records, manifest_rows, min_len, max_len, log):
    """Length filter applied to records AND their manifest rows in lockstep."""
    before = len(records)
    keep_ids = {r.id for r in records if min_len <= len(r.seq) <= max_len}
    f_recs = [r for r in records if r.id in keep_ids]
    f_man  = [m for m in manifest_rows if m["accession"] in keep_ids]
    if len(f_recs) < before:
        log.info(f"    Length filter [{min_len}-{max_len} bp]: {before} → {len(f_recs)}")
    return f_recs, f_man


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 2 — ALIGNMENT                                                      ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def detect_aligner() -> str:
    for tool in ("mafft", "muscle"):
        if shutil.which(tool):
            return tool
    if shutil.which("mafft.bat"):
        return "mafft"
    return "kmer"


def align_mafft(fasta_in, fasta_out, threads, log) -> bool:
    cmd = ["mafft", "--auto", "--thread", str(threads), "--quiet", fasta_in]
    try:
        with open(fasta_out, "w") as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.DEVNULL, check=True)
        log.info(f"  MAFFT alignment → {fasta_out}")
        return True
    except Exception as e:
        log.warning(f"  MAFFT failed: {e}")
        return False


def align_muscle(fasta_in, fasta_out, log) -> bool:
    for cmd in (
        ["muscle", "-align", fasta_in, "-output", fasta_out],   # v5
        ["muscle", "-in", fasta_in, "-out", fasta_out, "-quiet"], # v3
    ):
        try:
            subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)
            log.info(f"  MUSCLE alignment → {fasta_out}")
            return True
        except Exception:
            continue
    log.warning("  MUSCLE failed (tried v3 and v5 CLIs)")
    return False


def align_sequences(fasta_in, fasta_out, aligner, threads, log) -> bool:
    if aligner == "mafft":
        return align_mafft(fasta_in, fasta_out, threads, log)
    if aligner == "muscle":
        return align_muscle(fasta_in, fasta_out, log)
    shutil.copy(fasta_in, fasta_out)
    log.info(f"  k-mer mode: no alignment performed (copy → {fasta_out})")
    return True


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 3 — DISTANCE MATRIX                                                ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def encode_alignment(records) -> tuple:
    n = len(records)
    L = max(len(r.seq) for r in records)
    arr = np.zeros((n, L), dtype=np.uint8)
    code = {ord("A"): 1, ord("C"): 2, ord("G"): 3, ord("T"): 4, ord("U"): 4,
            ord("a"): 1, ord("c"): 2, ord("g"): 3, ord("t"): 4, ord("u"): 4}
    for i, r in enumerate(records):
        row = np.frombuffer(str(r.seq).encode("ascii", errors="replace"),
                            dtype=np.uint8)
        out = np.zeros(len(row), dtype=np.uint8)
        for k, v in code.items():
            out[row == k] = v
        arr[i, :len(out)] = out
    valid = arr > 0
    return arr, valid


def vectorised_p_distance(arr, valid, log):
    n = arr.shape[0]
    D = np.zeros((n, n), dtype=np.float64)
    log.info(f"  Vectorised p-distance for {n} sequences...")
    for i in tqdm(range(n), desc="  Distances", leave=False):
        ai = arr[i]; vi = valid[i]
        comp = valid[i + 1:] & vi
        diff = (arr[i + 1:] != ai) & comp
        comp_n = comp.sum(axis=1)
        diff_n = diff.sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            d = np.where(comp_n > 0, diff_n / comp_n, 1.0)
        D[i, i + 1:] = d
        D[i + 1:, i] = d
    return D


def jukes_cantor_correct(D):
    p = np.clip(D, 0.0, 0.7499)
    with np.errstate(invalid="ignore"):
        return -0.75 * np.log(1.0 - (4.0 / 3.0) * p)


def kmer_distance_matrix(records, k=6, log=None):
    def _km(s):
        s = s.upper().replace("-", "").replace("N", "")
        return set(s[i:i + k] for i in range(len(s) - k + 1))
    sets = [_km(str(r.seq)) for r in records]
    n = len(sets)
    D = np.zeros((n, n))
    it = range(n)
    if log:
        log.info(f"  k-mer Jaccard distance (n={n}, k={k})...")
        it = tqdm(it, desc="  Distances", leave=False)
    for i in it:
        for j in range(i + 1, n):
            si, sj = sets[i], sets[j]
            d = 1.0 - len(si & sj) / len(si | sj) if (si and sj) else 1.0
            D[i, j] = D[j, i] = d
    return D


def compute_distance_matrix(records, use_kmer, jc_correct, log):
    ids = [r.id for r in records]
    if use_kmer:
        D = kmer_distance_matrix(records, log=log)
        if jc_correct:
            log.warning("  --jc-correct ignored: not valid for k-mer Jaccard.")
    else:
        arr, valid = encode_alignment(records)
        D = vectorised_p_distance(arr, valid, log)
        if jc_correct:
            log.info("  Applying Jukes-Cantor correction...")
            D = jukes_cantor_correct(D)
            finite_max = np.nanmax(D[np.isfinite(D)]) if np.any(np.isfinite(D)) else 1.0
            D = np.where(np.isfinite(D), D, finite_max)
    return D, ids


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 4 — CLUSTERING & NOVEL VARIANT FLAGGING                            ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def _parse_tag(desc, key, default="Unknown"):
    m = re.search(rf"\[{re.escape(key)}=([^\]]+)\]", desc)
    return m.group(1).strip() if m else default


def _natural_st_key(st: str) -> tuple:
    """
    Sort key for subtype labels in natural numeric order:
      ST1, ST2, ..., ST9, ST10, ST11, ..., ST44,
      NMAST1, NMAST2, ...,
      novel_1, novel_2, ...,
      ST_unassigned, Unknown
    rather than alphabetical (which puts ST10 between ST1 and ST2).

    Group ordering: ST (0) → NMAST (1) → novel (2) → unassigned (8) → other (9).
    Within each group, numeric value drives the order.
    """
    if not isinstance(st, str):
        return (9, 0, "")
    if st.startswith("NMAST"):
        m = re.search(r"\d+", st)
        return (1, int(m.group()) if m else 0, st)
    if st.startswith("novel_") or st.startswith("novel"):
        m = re.search(r"\d+", st)
        return (2, int(m.group()) if m else 0, st)
    if st == "ST_unassigned":
        return (8, 0, st)
    if st.startswith("ST"):
        m = re.search(r"\d+", st)
        return (0, int(m.group()) if m else 0, st)
    return (9, 0, st)


def natsort_sts(sts):
    """Return a list of ST labels in natural numeric order."""
    return sorted(sts, key=_natural_st_key)


def cluster_sequences(D, ids, records, log,
                      thresh_cluster=0.05, thresh_subtype=0.025,
                      thresh_variant=0.01, novel_distance_multiplier=3.0):
    n = len(ids)
    if n < 2:
        log.warning("  <2 sequences — skipping clustering.")
        return pd.DataFrame({"accession": ids, "subtype": ["Unknown"]*n}), None

    Z = linkage(squareform(D, checks=False), method="average")
    cluster_cols = {
        "cluster_top":     fcluster(Z, t=thresh_cluster, criterion="distance"),
        "cluster_subtype": fcluster(Z, t=thresh_subtype, criterion="distance"),
        "cluster_variant": fcluster(Z, t=thresh_variant, criterion="distance"),
    }

    df = pd.DataFrame({
        "accession":       ids,
        "subtype":         [_parse_tag(r.description, "subtype")     for r in records],
        "host":            [_parse_tag(r.description, "host")        for r in records],
        "host_source":     [_parse_tag(r.description, "host_source") for r in records],
        "country":         [_parse_tag(r.description, "country")     for r in records],
        "collection_date": [_parse_tag(r.description, "date")        for r in records],
        "isolation_source":[_parse_tag(r.description, "isolation_source") for r in records],
        "seq_technology":  [_parse_tag(r.description, "seq_tech")    for r in records],
        "st_confidence":   [_parse_tag(r.description, "confidence")  for r in records],
        "seq_length":      [len(r.seq) for r in records],
        **cluster_cols,
    })

    def _status_for(s):
        if isinstance(s, str):
            if s.startswith("novel_"):
                return "novel"
            if s == "ST_unassigned":
                return "unassigned"
        return SUBTYPE_REGISTRY.get(s, {}).get("status", "unknown")
    df["st_status"]    = df["subtype"].map(_status_for)
    df["human_report"] = df["subtype"].map(lambda s: SUBTYPE_REGISTRY.get(s, {}).get("human", False))
    df["st_note"]      = df["subtype"].map(lambda s: SUBTYPE_REGISTRY.get(s, {}).get("note", ""))

    df["dist_from_medoid"]        = np.nan
    df["potential_novel_variant"] = False
    df["novel_variant_reason"]    = ""

    for st in df["subtype"].unique():
        if st == "Unknown":
            continue
        idx = df.index[df["subtype"] == st].to_numpy()
        if len(idx) < 3:
            continue

        sub_D = D[np.ix_(idx, idx)]
        medoid_pos = int(np.argmin(sub_D.sum(axis=1)))
        d_med = sub_D[medoid_pos]
        df.loc[idx, "dist_from_medoid"] = d_med

        triu = sub_D[np.triu_indices_from(sub_D, k=1)]
        median_intra = float(np.median(triu)) if len(triu) else 0.0
        scale = max(median_intra, 0.002)
        cutoff = max(novel_distance_multiplier * scale, 2 * thresh_variant)

        st_clusters = df.loc[idx, "cluster_subtype"]
        sizes = st_clusters.value_counts()
        dominant = sizes.idxmax()
        dominant_frac = sizes.max() / sizes.sum()

        for i, abs_idx in enumerate(idx):
            d = d_med[i]
            in_minority = (st_clusters.iloc[i] != dominant) and (dominant_frac > 0.6)
            reasons = []
            if d > cutoff:
                reasons.append(f"d_medoid={d:.4f} > {cutoff:.4f}")
            if in_minority:
                reasons.append(f"minority cluster ({st_clusters.iloc[i]})")
            if reasons:
                df.at[abs_idx, "potential_novel_variant"] = True
                df.at[abs_idx, "novel_variant_reason"] = "; ".join(reasons)

    novel = int(df["potential_novel_variant"].sum())
    log.info(f"  Clustering complete. {novel} potential novel-variant candidates.")
    return df, Z


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 5 — NEIGHBOR-JOINING TREE                                          ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def neighbor_joining(D, ids):
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    from Bio import Phylo
    dm = DistanceMatrix(ids, [[D[i][j] for j in range(i + 1)] for i in range(len(ids))])
    tree = DistanceTreeConstructor().nj(dm)
    h = StringIO()
    Phylo.write(tree, h, "newick")
    return h.getvalue()


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 6 — ORDINATION                                                     ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def pcoa(D):
    D2 = D ** 2
    row_mean = D2.mean(axis=1, keepdims=True)
    col_mean = D2.mean(axis=0, keepdims=True)
    grand = D2.mean()
    B = -0.5 * (D2 - row_mean - col_mean + grand)
    ev, evec = np.linalg.eigh(B)
    idx = np.argsort(ev)[::-1]
    ev = ev[idx]; evec = evec[:, idx]
    pos = ev > 0
    coords = evec[:, pos] * np.sqrt(ev[pos])
    var = ev[pos] / ev[pos].sum()
    return coords[:, :2], var[:2] if len(var) >= 2 else var


def compute_embedding(D, log, seed=42):
    if HAS_UMAP and D.shape[0] >= 10:
        log.info("  Computing UMAP embedding...")
        reducer = umap.UMAP(n_components=2, metric="precomputed",
                            n_neighbors=min(15, D.shape[0] - 1), random_state=seed)
        return reducer.fit_transform(D), "UMAP", None
    log.info("  Computing PCoA embedding...")
    coords, var = pcoa(D)
    return coords, "PCoA", var


def cliffs_delta(x, y):
    x = np.asarray(x); y = np.asarray(y)
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    cap = 5000
    if len(x) > cap:
        x = np.random.default_rng(0).choice(x, cap, replace=False)
    if len(y) > cap:
        y = np.random.default_rng(0).choice(y, cap, replace=False)
    gt = (x[:, None] > y[None, :]).sum()
    lt = (x[:, None] < y[None, :]).sum()
    return (gt - lt) / (len(x) * len(y))


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  STEP 7 — PLOTS & REPORTS                                                ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def plot_ordination(coords, cluster_df, method, var, out_dir, log):
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for st in natsort_sts(cluster_df["subtype"].unique()):
        mask = cluster_df["subtype"] == st
        idx = cluster_df.index[mask].to_numpy()
        status = SUBTYPE_REGISTRY.get(st, {}).get("status", "accepted")
        mk = STATUS_MARKERS.get(status, STATUS_MARKERS["accepted"])
        axes[0].scatter(coords[idx, 0], coords[idx, 1],
                        c=SUBTYPE_COLORS.get(st, "#888"), label=st, s=28,
                        alpha=mk["alpha"], marker=mk["marker"],
                        edgecolors=mk["edgecolor"])
    axes[0].set_title(f"{method} — by subtype", fontweight="bold")
    if method == "PCoA" and var is not None and len(var) >= 2:
        axes[0].set_xlabel(f"PCoA 1 ({var[0]*100:.1f}%)")
        axes[0].set_ylabel(f"PCoA 2 ({var[1]*100:.1f}%)")
    else:
        axes[0].set_xlabel(f"{method} Axis 1")
        axes[0].set_ylabel(f"{method} Axis 2")
    axes[0].legend(title="Subtype", bbox_to_anchor=(1.01, 1), loc="upper left",
                   fontsize=7, markerscale=1.5)
    axes[0].grid(alpha=0.2)

    novel = cluster_df["potential_novel_variant"].to_numpy()
    axes[1].scatter(coords[~novel, 0], coords[~novel, 1],
                    c="#AAAAAA", s=20, alpha=0.5, label="Typical", edgecolors="none")
    axes[1].scatter(coords[novel, 0], coords[novel, 1],
                    c="#E41A1C", s=50, alpha=0.9, label="Novel candidate",
                    edgecolors="black", linewidths=0.5, zorder=5)
    axes[1].set_title(f"{method} — novel-variant candidates", fontweight="bold")
    axes[1].set_xlabel(f"{method} Axis 1"); axes[1].set_ylabel(f"{method} Axis 2")
    axes[1].legend(); axes[1].grid(alpha=0.2)

    plt.tight_layout()
    p = os.path.join(out_dir, "plots", "ordination.png")
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    log.info(f"  ✔ {p}")


def plot_distance_heatmap(D, cluster_df, out_dir, log):
    order = cluster_df.sort_values(["subtype", "cluster_subtype"]).index.to_numpy()
    D_sorted = D[np.ix_(order, order)]
    fig, ax = plt.subplots(figsize=(12, 10))
    vmax = np.quantile(D_sorted[np.triu_indices_from(D_sorted, k=1)], 0.95)
    im = ax.imshow(D_sorted, cmap="viridis_r", aspect="auto", vmin=0, vmax=vmax)
    plt.colorbar(im, ax=ax, label="pairwise distance", fraction=0.03)

    sts = cluster_df.loc[order, "subtype"].values
    boundaries = [0]
    for i in range(1, len(sts)):
        if sts[i] != sts[i - 1]:
            boundaries.append(i)
    boundaries.append(len(sts))
    for b in boundaries[1:-1]:
        ax.axhline(b - 0.5, color="white", lw=1.5)
        ax.axvline(b - 0.5, color="white", lw=1.5)
    ticks = [(boundaries[i] + boundaries[i + 1]) // 2
             for i in range(len(boundaries) - 1)]
    labels = [sts[t] for t in ticks]
    ax.set_xticks(ticks); ax.set_xticklabels(labels, rotation=45, fontsize=8)
    ax.set_yticks(ticks); ax.set_yticklabels(labels, fontsize=8)
    ax.set_title(f"Pairwise distance heatmap (sorted by ST; scale capped at 95th pct = {vmax:.3f})",
                 fontweight="bold")
    plt.tight_layout()
    p = os.path.join(out_dir, "plots", "distance_heatmap.png")
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    log.info(f"  ✔ {p}")


def _is_known_st(label: str) -> bool:
    """Classify a subtype label: True for accepted STs, False for novel_N / unassigned."""
    if not isinstance(label, str):
        return False
    if label.startswith("novel_") or label == "ST_unassigned":
        return False
    return label.startswith("ST") or label.startswith("NMAST")


def plot_intra_inter(D, cluster_df, out_dir, log, thresh_subtype=0.025):
    """
    Two figures:
      intra_inter_boxplot.png         — overall intra vs inter + per-ST (known STs only)
      intra_inter_novel.png           — per-cluster intra distance for novel_N + ST_unassigned
    Separating the panels prevents the y-axis scale on the main figure from being
    blown out by novel clusters that contain divergent records.
    """
    sts = cluster_df["subtype"].values
    n = len(sts)
    iu, ju = np.triu_indices(n, k=1)
    same = sts[iu] == sts[ju]
    not_unk = (sts[iu] != "Unknown") & (sts[ju] != "Unknown")
    d = D[iu, ju]

    # For the main panel: restrict to known STs so a divergent novel cluster
    # doesn't compress the y-axis
    known_mask_i = np.array([_is_known_st(s) for s in sts])
    both_known = known_mask_i[iu] & known_mask_i[ju]

    intra = d[same & not_unk & both_known].tolist()
    inter = d[~same & not_unk & both_known].tolist()

    st_intra_known = defaultdict(list)
    st_intra_novel = defaultdict(list)
    sel = same & not_unk
    for i_, dv in zip(iu[sel], d[sel]):
        s = sts[i_]
        if _is_known_st(s):
            st_intra_known[s].append(float(dv))
        else:
            st_intra_novel[s].append(float(dv))

    # ── PANEL 1 — Known STs only ─────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    bp = axes[0].boxplot([intra, inter], patch_artist=True, notch=True,
                         tick_labels=["Intra-subtype", "Inter-subtype"])
    bp["boxes"][0].set_facecolor("#4DAF4A")
    bp["boxes"][1].set_facecolor("#E41A1C")
    for e in ["whiskers", "caps", "medians", "fliers"]:
        for it in bp[e]:
            it.set_color("black")
    axes[0].set_ylabel("Distance"); axes[0].grid(axis="y", alpha=0.3)
    axes[0].set_title("Intra- vs inter-subtype distances\n(known STs only)",
                      fontweight="bold")

    if intra and inter:
        delta = cliffs_delta(intra, inter)
        axes[0].text(
            1.5, max(max(intra), max(inter)) * 0.95,
            f"median intra = {np.median(intra):.4f}\n"
            f"median inter = {np.median(inter):.4f}\n"
            f"Cliff's δ = {delta:.3f}",
            ha="center", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="gray"),
        )

    sorted_sts = natsort_sts(st_intra_known.keys())
    data = [st_intra_known[st] for st in sorted_sts if st_intra_known[st]]
    labels = [st for st in sorted_sts if st_intra_known[st]]
    colors = [SUBTYPE_COLORS.get(st, "#888") for st in labels]
    bp2 = axes[1].boxplot(data, patch_artist=True, tick_labels=labels)
    for patch, c in zip(bp2["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.8)
    axes[1].set_ylabel("Intra-subtype distance")
    axes[1].set_title("Intra-subtype distance per known subtype", fontweight="bold")
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=45)
    axes[1].axhline(thresh_subtype, color="red", ls="--", lw=1.5,
                    label=f"Subtype threshold ({thresh_subtype})")
    axes[1].legend(fontsize=8); axes[1].grid(axis="y", alpha=0.3)

    # Cap y-axis at the 99th percentile of all known-ST intra distances if
    # there are obvious outliers above ~0.3 — these are almost always
    # contaminating records (RefSeq predicted, mis-assigned ST etc.) that
    # should be investigated, not used to set the plot scale
    all_intra = [v for vs in st_intra_known.values() for v in vs]
    if all_intra:
        p99 = np.quantile(all_intra, 0.99)
        if p99 < 0.2:
            axes[1].set_ylim(-0.005, max(0.15, p99 * 1.5))

    plt.tight_layout()
    p = os.path.join(out_dir, "plots", "intra_inter_boxplot.png")
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    log.info(f"  ✔ {p}")

    # ── PANEL 2 — novel_N and ST_unassigned, separately ─────────────────
    if st_intra_novel:
        fig2, ax2 = plt.subplots(figsize=(max(8, 0.5 * len(st_intra_novel)), 6))
        sorted_novel = natsort_sts(st_intra_novel.keys())
        data2 = [st_intra_novel[st] for st in sorted_novel if st_intra_novel[st]]
        labels2 = [st for st in sorted_novel if st_intra_novel[st]]
        colors2 = [SUBTYPE_COLORS.get(st, "#888") for st in labels2]
        bp3 = ax2.boxplot(data2, patch_artist=True, tick_labels=labels2)
        for patch, c in zip(bp3["boxes"], colors2):
            patch.set_facecolor(c); patch.set_alpha(0.8)
        ax2.set_ylabel("Intra-cluster distance")
        ax2.set_title("Intra-cluster distance — novel candidates and unassigned\n"
                      "(separated from known STs; note the larger y-axis range)",
                      fontweight="bold")
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha="right")
        ax2.axhline(thresh_subtype, color="red", ls="--", lw=1.5,
                    label=f"Known-ST threshold ({thresh_subtype})")
        ax2.axhline(0.20, color="orange", ls="--", lw=1.5,
                    label="Novel-cluster cutoff (0.20)")
        ax2.legend(fontsize=8, loc="upper left")
        ax2.grid(axis="y", alpha=0.3)
        plt.tight_layout()
        p2 = os.path.join(out_dir, "plots", "intra_inter_novel.png")
        plt.savefig(p2, dpi=150, bbox_inches="tight"); plt.close()
        log.info(f"  ✔ {p2}")


def plot_novel_clusters(cluster_df, out_dir, log):
    """
    Dedicated four-panel figure for the novel_N provisional clusters,
    showing for each cluster: size, host distribution, country distribution,
    sequencing-technology distribution, and collection-year coverage.

    Helps the reader judge whether a novel_N cluster is "a coherent
    single-study artefact" (one host + one country + one tech + one year
    suggests a single deposit batch) or "a genuinely distinct lineage"
    (spread across multiple hosts/countries/years).
    """
    df = cluster_df.copy()
    df = df[df["subtype"].str.startswith("novel_", na=False)]
    if df.empty:
        log.info("  No novel clusters — skipping novel-cluster facet plot.")
        return

    df["country_clean"] = df["country"].astype(str).map(_normalise_country)
    df["year"] = df["collection_date"].astype(str).map(_extract_year)

    novel_labels = natsort_sts(df["subtype"].unique())

    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28)
    ax_size  = fig.add_subplot(gs[0, 0])
    ax_host  = fig.add_subplot(gs[0, 1])
    ax_geo   = fig.add_subplot(gs[1, 0])
    ax_tech  = fig.add_subplot(gs[1, 1])

    # ── PANEL A — Cluster sizes ─────────────────────────────────────────
    sizes = df["subtype"].value_counts().reindex(novel_labels).fillna(0).astype(int)
    bars = ax_size.bar(range(len(sizes)), sizes.values,
                       color=[SUBTYPE_COLORS.get(l, "#888") for l in sizes.index],
                       edgecolor="white", linewidth=0.8)
    ax_size.set_xticks(range(len(sizes)))
    ax_size.set_xticklabels(sizes.index, rotation=45, ha="right", fontsize=9)
    ax_size.set_ylabel("records in cluster")
    ax_size.set_title("Novel cluster sizes", fontweight="bold")
    for i, v in enumerate(sizes.values):
        ax_size.text(i, v + 0.1, str(v), ha="center", fontsize=9)
    ax_size.grid(axis="y", alpha=0.3)

    # ── PANEL B — Host × novel_N heatmap ────────────────────────────────
    host_df = df[df["host"].astype(str) != "Unknown"]
    if not host_df.empty:
        host_pivot = (host_df.groupby(["host", "subtype"]).size()
                      .unstack(fill_value=0))
        host_pivot = host_pivot.reindex(columns=novel_labels, fill_value=0)
        # Sort hosts by total records descending
        host_pivot = host_pivot.loc[
            host_pivot.sum(axis=1).sort_values(ascending=False).index]
        # Limit to top 12 hosts so the heatmap is readable
        host_pivot = host_pivot.head(12)
        if not host_pivot.empty and host_pivot.values.sum() > 0:
            sns.heatmap(host_pivot, annot=True, fmt="d", cmap="YlGnBu",
                        linewidths=0.4, ax=ax_host, cbar_kws={"label": "count"},
                        annot_kws={"size": 9})
            ax_host.set_title("Host × novel cluster", fontweight="bold")
            ax_host.set_xlabel(""); ax_host.set_ylabel("Host")
            ax_host.tick_params(axis="x", rotation=45, labelsize=9)
            ax_host.tick_params(axis="y", labelsize=9)
    if df[df["host"].astype(str) != "Unknown"].empty:
        ax_host.text(0.5, 0.5, "no host data for novel clusters",
                     ha="center", va="center", fontsize=11, color="#888",
                     transform=ax_host.transAxes)
        ax_host.axis("off")

    # ── PANEL C — Country × novel_N heatmap ─────────────────────────────
    geo_df = df[df["country_clean"] != "Unknown"]
    if not geo_df.empty:
        geo_pivot = (geo_df.groupby(["country_clean", "subtype"]).size()
                     .unstack(fill_value=0))
        geo_pivot = geo_pivot.reindex(columns=novel_labels, fill_value=0)
        geo_pivot = geo_pivot.loc[
            geo_pivot.sum(axis=1).sort_values(ascending=False).index]
        geo_pivot = geo_pivot.head(12)
        if not geo_pivot.empty and geo_pivot.values.sum() > 0:
            sns.heatmap(geo_pivot, annot=True, fmt="d", cmap="Reds",
                        linewidths=0.4, ax=ax_geo, cbar_kws={"label": "count"},
                        annot_kws={"size": 9})
            ax_geo.set_title("Country × novel cluster", fontweight="bold")
            ax_geo.set_xlabel(""); ax_geo.set_ylabel("Country")
            ax_geo.tick_params(axis="x", rotation=45, labelsize=9)
            ax_geo.tick_params(axis="y", labelsize=9)
    if geo_df.empty:
        ax_geo.text(0.5, 0.5, "no country data for novel clusters",
                    ha="center", va="center", fontsize=11, color="#888",
                    transform=ax_geo.transAxes)
        ax_geo.axis("off")

    # ── PANEL D — Sequencing technology × novel_N stacked bar ──────────
    if "seq_technology" in df.columns:
        tech_df = df.copy()
        tech_df["seq_technology"] = tech_df["seq_technology"].fillna("Unknown")
        tech_pivot = (tech_df.groupby(["subtype", "seq_technology"]).size()
                      .unstack(fill_value=0))
        tech_pivot = tech_pivot.reindex(index=novel_labels, fill_value=0)
        # Order tech columns by their global frequency (most common at bottom)
        tech_order = [t for t, _ in TECH_PATTERNS] + ["Unknown"]
        tech_pivot = tech_pivot.reindex(
            columns=[t for t in tech_order if t in tech_pivot.columns],
            fill_value=0)
        bar_colors = [TECH_COLORS.get(t, "#AAAAAA") for t in tech_pivot.columns]
        tech_pivot.plot(kind="bar", stacked=True, ax=ax_tech,
                        color=bar_colors, edgecolor="white", linewidth=0.5,
                        width=0.85)
        ax_tech.set_xlabel(""); ax_tech.set_ylabel("records")
        ax_tech.set_title("Sequencing technology × novel cluster",
                          fontweight="bold")
        ax_tech.tick_params(axis="x", rotation=45, labelsize=9)
        ax_tech.legend(title="technology", bbox_to_anchor=(1.01, 1),
                       loc="upper left", fontsize=8)
        ax_tech.grid(axis="y", alpha=0.3)

    # Suptitle with summary
    n_novel_records = len(df)
    n_clusters = df["subtype"].nunique()
    n_singletons = ((df["subtype"] == "ST_unassigned").sum()
                    if "ST_unassigned" in cluster_df["subtype"].values else 0)
    fig.suptitle(
        f"Novel-candidate clusters: {n_clusters} groups, "
        f"{n_novel_records} records "
        f"(provisional labels; not formal subtype proposals)",
        fontsize=13, fontweight="bold", y=0.995,
    )
    p = os.path.join(out_dir, "plots", "novel_clusters.png")
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    log.info(f"  ✔ {p}")


def plot_technology_breakdown(cluster_df, out_dir, log):
    if "seq_technology" not in cluster_df.columns:
        return
    df = cluster_df.copy()
    df["seq_technology"] = df["seq_technology"].fillna("Unknown")
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    tech_counts = df["seq_technology"].value_counts()
    colors = [TECH_COLORS.get(t, "#AAAAAA") for t in tech_counts.index]
    wedges, _, autotexts = axes[0].pie(
        tech_counts.values, labels=tech_counts.index, colors=colors,
        autopct=lambda p: f"{p:.1f}%" if p > 2 else "",
        startangle=140, wedgeprops=dict(edgecolor="white", linewidth=1.5),
    )
    for at in autotexts:
        at.set_fontsize(8)
    axes[0].set_title("Sequencing technology distribution", fontweight="bold")

    # Technology per subtype — known STs only (novel clusters get their own plot)
    df_known_st_tech = df[df["subtype"].map(_is_known_st)]
    pivot = df_known_st_tech.groupby(["subtype", "seq_technology"]).size().unstack(fill_value=0)
    pivot = pivot.loc[natsort_sts(pivot.index)]
    tech_order = [t for t, _ in TECH_PATTERNS] + ["Unknown"]
    pivot = pivot.reindex(columns=[c for c in tech_order if c in pivot.columns],
                          fill_value=0)
    bar_colors = [TECH_COLORS.get(t, "#AAAAAA") for t in pivot.columns]
    pivot.plot(kind="bar", stacked=True, ax=axes[1],
               color=bar_colors, edgecolor="white", linewidth=0.5)
    axes[1].set_title("Technology per known subtype", fontweight="bold")
    axes[1].set_xlabel("Subtype"); axes[1].set_ylabel("Sequence count")
    axes[1].tick_params(axis="x", rotation=45)
    axes[1].legend(title="Technology", bbox_to_anchor=(1.01, 1), loc="upper left",
                   fontsize=8)
    axes[1].grid(axis="y", alpha=0.3)

    hm = pivot.T
    hm = hm.loc[(hm > 0).any(axis=1)]
    if not hm.empty:
        sns.heatmap(hm, annot=True, fmt="d", cmap="Blues",
                    linewidths=0.5, linecolor="lightgray", ax=axes[2],
                    cbar_kws={"label": "count"})
        axes[2].set_title("Technology × subtype", fontweight="bold")
        axes[2].set_xlabel("Subtype"); axes[2].set_ylabel("Technology")
        axes[2].tick_params(axis="x", rotation=45)

    plt.tight_layout()
    p = os.path.join(out_dir, "plots", "technology_breakdown.png")
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    log.info(f"  ✔ {p}")


def plot_host_distribution(cluster_df, out_dir, log):
    """Host frequency by known subtype + overall donut. Novel clusters are
    plotted separately in plot_novel_clusters."""
    if "host" not in cluster_df.columns:
        return
    df = cluster_df.copy()
    df_known = df[df["host"] != "Unknown"]
    if df_known.empty:
        log.info("  No known hosts — skipping host plot.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Overall donut — all records with known host
    host_counts = df_known["host"].value_counts().head(15)
    other_n = len(df_known) - host_counts.sum()
    if other_n > 0:
        host_counts["Other"] = other_n
    axes[0].pie(host_counts.values, labels=host_counts.index,
                autopct=lambda p: f"{p:.1f}%" if p > 3 else "",
                startangle=90, wedgeprops=dict(width=0.4, edgecolor="white"))
    axes[0].set_title("Top hosts (overall, excl. Unknown)", fontweight="bold")

    # Host x subtype heatmap — known STs only (novel/unassigned have a dedicated plot)
    df_known_st = df_known[df_known["subtype"].map(_is_known_st)]
    top_hosts = df_known_st["host"].value_counts().head(12).index
    sub_df = df_known_st[df_known_st["host"].isin(top_hosts)]
    pivot = sub_df.groupby(["host", "subtype"]).size().unstack(fill_value=0)
    if not pivot.empty:
        # Reorder columns by natural ST order (ST1, ST2, ..., ST10, ..., NMAST1, ...)
        # and rows by descending row sum for visual emphasis
        pivot = pivot.reindex(columns=natsort_sts(pivot.columns))
        pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).index]
        sns.heatmap(pivot, annot=True, fmt="d", cmap="YlOrRd",
                    linewidths=0.5, ax=axes[1],
                    cbar_kws={"label": "count"})
        axes[1].set_title("Host × subtype (known STs, top 12 hosts)",
                          fontweight="bold")
        axes[1].set_xlabel("Subtype"); axes[1].set_ylabel("Host")

    plt.tight_layout()
    p = os.path.join(out_dir, "plots", "host_distribution.png")
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    log.info(f"  ✔ {p}")


def _normalise_country(s: str) -> str:
    """
    Country values from GenBank /geo_loc_name look like 'Spain', 'Spain:Madrid',
    'France:Cote d'Azur, Antibes', 'Iraq: Erbil'. For aggregation, strip
    the trailing region/locality after the first colon and trim whitespace.
    """
    if not isinstance(s, str) or not s or s == "Unknown":
        return "Unknown"
    head = s.split(":", 1)[0].strip()
    # Some submitters use a comma instead of a colon
    head = head.split(",", 1)[0].strip()
    return head if head else "Unknown"


def _extract_year(s: str):
    """
    Parse a /collection_date string into a 4-digit year or None.
    Handles 'YYYY', 'Mon-YYYY', 'DD-Mon-YYYY', 'YYYY-MM-DD', 'YYYY/YYYY' (range
    — take latest), and similar variants.
    """
    if not isinstance(s, str) or not s or s == "Unknown":
        return None
    # Find all 4-digit years (1900-2099) in the string and take the largest
    # (handles ranges like '1998/2021' by preferring the later year)
    years = re.findall(r"(19[0-9]{2}|20[0-9]{2})", s)
    if not years:
        return None
    return int(max(years))


def plot_geographic_temporal(cluster_df, out_dir, log):
    """
    Four-panel plot of country & temporal distribution:
      (a) Top countries overall (horizontal bar, log scale)
      (b) Country × subtype heatmap (top 12 countries, all STs incl. novel/unassigned)
      (c) Records per year (overall + stacked by ST status)
      (d) Country × year heatmap (when did each country contribute records?)
    """
    if "country" not in cluster_df.columns or "collection_date" not in cluster_df.columns:
        log.info("  Skipping geographic_temporal plot (missing columns).")
        return

    df = cluster_df.copy()
    df["country_clean"] = df["country"].astype(str).map(_normalise_country)
    df["year"] = df["collection_date"].astype(str).map(_extract_year)
    df_known_country = df[df["country_clean"] != "Unknown"]
    df_known_year    = df[df["year"].notna()]

    n_unknown_country = (df["country_clean"] == "Unknown").sum()
    n_unknown_year    = df["year"].isna().sum()
    log.info(f"  Country: {len(df) - n_unknown_country}/{len(df)} known; "
             f"Year: {len(df) - n_unknown_year}/{len(df)} known")

    fig = plt.figure(figsize=(20, 14))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.35, wspace=0.25)
    ax_country = fig.add_subplot(gs[0, 0])
    ax_cmap    = fig.add_subplot(gs[0, 1])
    ax_year    = fig.add_subplot(gs[1, 0])
    ax_cy      = fig.add_subplot(gs[1, 1])

    # (a) Top countries — horizontal bar
    country_counts = df_known_country["country_clean"].value_counts().head(15)
    if not country_counts.empty:
        bars = ax_country.barh(
            range(len(country_counts)), country_counts.values,
            color="#4DAF4A", edgecolor="white",
        )
        ax_country.set_yticks(range(len(country_counts)))
        ax_country.set_yticklabels(country_counts.index, fontsize=9)
        ax_country.invert_yaxis()
        ax_country.set_xlabel("records")
        ax_country.set_title("Top countries (excl. Unknown)", fontweight="bold")
        for i, v in enumerate(country_counts.values):
            ax_country.text(v + max(country_counts.values) * 0.01, i,
                            str(v), va="center", fontsize=9)
        ax_country.grid(axis="x", alpha=0.3)

    # (b) Country × subtype heatmap
    if not df_known_country.empty:
        top_countries = df_known_country["country_clean"].value_counts().head(12).index
        sub = df_known_country[df_known_country["country_clean"].isin(top_countries)]
        pivot = sub.groupby(["country_clean", "subtype"]).size().unstack(fill_value=0)
        if not pivot.empty:
            pivot = pivot.reindex(columns=natsort_sts(pivot.columns))
            pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).index]
            sns.heatmap(pivot, annot=True, fmt="d", cmap="YlGnBu",
                        linewidths=0.5, ax=ax_cmap,
                        cbar_kws={"label": "count"})
            ax_cmap.set_title("Country × subtype (top 12 countries)", fontweight="bold")
            ax_cmap.set_xlabel("Subtype"); ax_cmap.set_ylabel("Country")
            ax_cmap.tick_params(axis="x", rotation=45, labelsize=9)
            ax_cmap.tick_params(axis="y", labelsize=9)

    # (c) Records per year, stacked by ST status group
    if not df_known_year.empty:
        df_known_year = df_known_year.copy()
        df_known_year["year"] = df_known_year["year"].astype(int)
        # Bucket STs into broad groups for the stack
        def _group_st(s):
            if isinstance(s, str):
                if s.startswith("novel_"):    return "novel candidates"
                if s == "ST_unassigned":      return "ST_unassigned"
                if s in ("ST1", "ST2", "ST3"): return "ST1-ST3 (dominant human)"
                if s in ("ST4", "ST5", "ST6", "ST7", "ST8", "ST9"):
                    return "ST4-ST9 (frequent)"
            return "other STs"
        df_known_year["st_group"] = df_known_year["subtype"].map(_group_st)
        year_pivot = (df_known_year.groupby(["year", "st_group"]).size()
                      .unstack(fill_value=0))
        # Preserve group order
        group_order = ["ST1-ST3 (dominant human)", "ST4-ST9 (frequent)",
                       "other STs", "novel candidates", "ST_unassigned"]
        group_order = [g for g in group_order if g in year_pivot.columns]
        year_pivot = year_pivot[group_order]
        # Fill any missing years between min and max so the x-axis is regular
        full_range = list(range(int(year_pivot.index.min()),
                                 int(year_pivot.index.max()) + 1))
        year_pivot = year_pivot.reindex(full_range, fill_value=0)
        colors = {"ST1-ST3 (dominant human)": "#E41A1C",
                  "ST4-ST9 (frequent)":       "#377EB8",
                  "other STs":                "#666666",
                  "novel candidates":         "#FF1493",
                  "ST_unassigned":            "#BBBBBB"}
        bottom = np.zeros(len(year_pivot))
        for grp in group_order:
            vals = year_pivot[grp].to_numpy()
            ax_year.bar(year_pivot.index, vals, bottom=bottom,
                        label=grp, color=colors[grp],
                        edgecolor="white", linewidth=0.4, width=0.85)
            bottom += vals
        ax_year.set_xlabel("Collection year")
        ax_year.set_ylabel("records")
        ax_year.set_title("Records per collection year", fontweight="bold")
        ax_year.legend(fontsize=8, loc="upper left")
        ax_year.grid(axis="y", alpha=0.3)
        if len(full_range) > 15:
            # Thin x-axis labels
            ax_year.set_xticks([y for y in full_range if y % 2 == 0])
        ax_year.tick_params(axis="x", rotation=45, labelsize=9)

    # (d) Country × year heatmap
    cy_df = df[(df["country_clean"] != "Unknown") & (df["year"].notna())].copy()
    if not cy_df.empty:
        cy_df["year"] = cy_df["year"].astype(int)
        top_c = cy_df["country_clean"].value_counts().head(10).index
        cy_top = cy_df[cy_df["country_clean"].isin(top_c)]
        cy_pivot = (cy_top.groupby(["country_clean", "year"]).size()
                    .unstack(fill_value=0))
        # Reindex year range
        if not cy_pivot.empty:
            yrs = list(range(int(cy_pivot.columns.min()),
                              int(cy_pivot.columns.max()) + 1))
            cy_pivot = cy_pivot.reindex(columns=yrs, fill_value=0)
            # Country order by total descending
            cy_pivot = cy_pivot.loc[
                cy_pivot.sum(axis=1).sort_values(ascending=False).index]
            sns.heatmap(cy_pivot, annot=True, fmt="d", cmap="Reds",
                        linewidths=0.3, ax=ax_cy,
                        cbar_kws={"label": "count"},
                        annot_kws={"size": 8})
            ax_cy.set_title("Country × collection year (top 10 countries)",
                            fontweight="bold")
            ax_cy.set_xlabel("Year"); ax_cy.set_ylabel("Country")
            ax_cy.tick_params(axis="x", rotation=45, labelsize=9)
            ax_cy.tick_params(axis="y", labelsize=9)

    # Annotate the unknown-metadata rates as a small note
    note = (f"records with unknown country: {n_unknown_country}/{len(df)} "
            f"({100*n_unknown_country/len(df):.1f}%)  ·  "
            f"records with unknown year: {n_unknown_year}/{len(df)} "
            f"({100*n_unknown_year/len(df):.1f}%)")
    fig.suptitle("Geographic and temporal distribution", fontsize=14,
                 fontweight="bold", y=0.995)
    fig.text(0.5, 0.005, note, ha="center", fontsize=9,
             style="italic", color="#444")

    p = os.path.join(out_dir, "plots", "geographic_temporal.png")
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    log.info(f"  ✔ {p}")


def plot_acquisition_breakdown(manifest_df, out_dir, log):
    """How was each record acquired? (pass1_narrow, pass2_reference, pass2_placement)"""
    if "acquisition_pass" not in manifest_df.columns:
        return
    if "st_confidence" not in manifest_df.columns:
        return

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    counts = manifest_df["acquisition_pass"].value_counts()
    axes[0].bar(range(len(counts)), counts.values,
                color=["#4DAF4A", "#377EB8", "#FF7F00"][:len(counts)])
    axes[0].set_xticks(range(len(counts)))
    axes[0].set_xticklabels(counts.index, rotation=15)
    axes[0].set_ylabel("Records"); axes[0].grid(axis="y", alpha=0.3)
    axes[0].set_title("Acquisition pass", fontweight="bold")
    for i, v in enumerate(counts.values):
        axes[0].text(i, v + 1, str(v), ha="center", fontsize=10)

    # Confidence breakdown
    conf_pivot = (manifest_df.groupby(["acquisition_pass", "st_confidence"])
                  .size().unstack(fill_value=0))
    conf_pivot.plot(kind="barh", stacked=True, ax=axes[1],
                    colormap="viridis", edgecolor="white")
    axes[1].set_xlabel("Records"); axes[1].set_ylabel("")
    axes[1].set_title("Confidence breakdown by acquisition pass",
                      fontweight="bold")
    axes[1].legend(title="Confidence", bbox_to_anchor=(1.01, 1), loc="upper left",
                   fontsize=8)
    axes[1].grid(axis="x", alpha=0.3)

    plt.tight_layout()
    p = os.path.join(out_dir, "plots", "acquisition_breakdown.png")
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    log.info(f"  ✔ {p}")


def build_interactive_report(coords, cluster_df, D, method, out_dir, log):
    plot_df = cluster_df.copy().reset_index(drop=True)
    plot_df["x"] = coords[:, 0]; plot_df["y"] = coords[:, 1]
    plot_df["hover"] = (
        "Accession: "  + plot_df["accession"] + "<br>" +
        "Subtype: "    + plot_df["subtype"]   + "<br>" +
        "Confidence: " + plot_df["st_confidence"].astype(str) + "<br>" +
        "Status: "     + plot_df["st_status"] + "<br>" +
        "Host: "       + plot_df["host"].astype(str) +
        " (" + plot_df["host_source"].astype(str) + ")<br>" +
        "Country: "    + plot_df["country"].astype(str) + "<br>" +
        "Tech: "       + plot_df["seq_technology"].astype(str) + "<br>" +
        "Cluster: "    + plot_df["cluster_variant"].astype(str) + "<br>" +
        "Novel: "      + plot_df["novel_variant_reason"].astype(str)
    )

    fig = make_subplots(
        rows=3, cols=2,
        subplot_titles=(
            f"{method} — by subtype", f"{method} — novel candidates",
            "Intra-subtype distances", "Sequence counts per subtype",
            "Sequencing technology", "Technology × subtype",
        ),
        specs=[[{"type": "xy"}, {"type": "xy"}],
               [{"type": "xy"}, {"type": "xy"}],
               [{"type": "domain"}, {"type": "xy"}]],
        vertical_spacing=0.10, horizontal_spacing=0.08,
    )

    for st in natsort_sts(plot_df["subtype"].unique()):
        sub = plot_df[plot_df["subtype"] == st]
        fig.add_trace(go.Scatter(
            x=sub["x"], y=sub["y"], mode="markers", name=st,
            marker=dict(color=SUBTYPE_COLORS.get(st, "#888"), size=6, opacity=0.8),
            text=sub["hover"], hoverinfo="text", legendgroup=st,
        ), row=1, col=1)

    typical = plot_df[~plot_df["potential_novel_variant"]]
    novel   = plot_df[plot_df["potential_novel_variant"]]
    fig.add_trace(go.Scatter(
        x=typical["x"], y=typical["y"], mode="markers", name="Typical",
        marker=dict(color="#AAAAAA", size=5, opacity=0.5),
        text=typical["hover"], hoverinfo="text", legendgroup="v",
    ), row=1, col=2)
    fig.add_trace(go.Scatter(
        x=novel["x"], y=novel["y"], mode="markers", name="⚠ Novel candidate",
        marker=dict(color="#E41A1C", size=10, opacity=0.95,
                    line=dict(color="black", width=1)),
        text=novel["hover"], hoverinfo="text", legendgroup="v",
    ), row=1, col=2)

    sts_arr = cluster_df["subtype"].values
    n = len(sts_arr)
    iu, ju = np.triu_indices(n, k=1)
    same = sts_arr[iu] == sts_arr[ju]
    not_unk = (sts_arr[iu] != "Unknown") & (sts_arr[ju] != "Unknown")
    st_intra = defaultdict(list)
    sel = same & not_unk
    for i_, dij in zip(iu[sel], D[iu[sel], ju[sel]]):
        st_intra[sts_arr[i_]].append(float(dij))
    for st in natsort_sts(st_intra.keys()):
        if st_intra[st]:
            fig.add_trace(go.Box(
                y=st_intra[st], name=st,
                marker_color=SUBTYPE_COLORS.get(st, "#888"),
                showlegend=False, boxpoints="outliers",
            ), row=2, col=1)

    counts = plot_df["subtype"].value_counts()
    counts = counts.reindex(natsort_sts(counts.index))
    fig.add_trace(go.Bar(
        x=counts.index.tolist(), y=counts.values.tolist(),
        marker_color=[SUBTYPE_COLORS.get(st, "#888") for st in counts.index],
        showlegend=False,
        hovertemplate="<b>%{x}</b><br>n=%{y}<extra></extra>",
    ), row=2, col=2)

    if "seq_technology" in cluster_df.columns:
        tech_counts = cluster_df["seq_technology"].fillna("Unknown").value_counts()
        fig.add_trace(go.Pie(
            labels=tech_counts.index.tolist(), values=tech_counts.values.tolist(),
            marker_colors=[TECH_COLORS.get(t, "#AAAAAA") for t in tech_counts.index],
            hole=0.35,
        ), row=3, col=1)

        tp = cluster_df.groupby(["subtype", "seq_technology"]).size().reset_index(name="n")
        for tech in tp["seq_technology"].unique():
            sub = tp[tp["seq_technology"] == tech].sort_values("subtype")
            fig.add_trace(go.Bar(
                x=sub["subtype"], y=sub["n"], name=tech,
                marker_color=TECH_COLORS.get(tech, "#AAAAAA"),
                legendgroup=f"tech_{tech}",
            ), row=3, col=2)

    fig.update_layout(
        barmode="stack",
        title_text="<b>Blastocystis SSU rRNA — v5 broad-first analysis</b>",
        title_x=0.5, height=1300,
        font=dict(family="Arial", size=11), paper_bgcolor="white",
    )

    p = os.path.join(out_dir, "plots", "interactive_report.html")
    fig.write_html(p)
    log.info(f"  ✔ {p}")


def write_summary_table(D, cluster_df, out_dir, log):
    sts = cluster_df["subtype"].values
    rows = []
    for st in natsort_sts(set(sts) - {"Unknown"}):
        idx = np.where(sts == st)[0]
        if len(idx) < 2:
            continue
        sub_D = D[np.ix_(idx, idx)]
        triu = sub_D[np.triu_indices_from(sub_D, k=1)]
        st_df = cluster_df.iloc[idx]

        tech_summary = "N/A"
        if "seq_technology" in st_df.columns:
            tc = st_df["seq_technology"].value_counts().to_dict()
            tech_summary = "; ".join(f"{t}:{c}" for t, c in
                                     sorted(tc.items(), key=lambda x: -x[1]))
        conf_summary = "N/A"
        if "st_confidence" in st_df.columns:
            cc = st_df["st_confidence"].value_counts().to_dict()
            conf_summary = "; ".join(f"{k}:{v}" for k, v in
                                     sorted(cc.items(), key=lambda x: -x[1]))
        host_summary = "N/A"
        if "host" in st_df.columns:
            hc = st_df["host"].value_counts().head(3).to_dict()
            host_summary = "; ".join(f"{h}:{c}" for h, c in hc.items())

        reg = SUBTYPE_REGISTRY.get(st, {})
        rows.append({
            "subtype": st,
            "status": reg.get("status", "unknown"),
            "reported_in_humans": reg.get("human", False),
            "registry_note": reg.get("note", ""),
            "n_sequences": len(idx),
            "mean_intra":   round(float(np.mean(triu)), 5)   if len(triu) else None,
            "median_intra": round(float(np.median(triu)), 5) if len(triu) else None,
            "max_intra":    round(float(np.max(triu)), 5)    if len(triu) else None,
            "n_variant_clusters": int(st_df["cluster_variant"].nunique()),
            "n_novel_candidates": int(st_df["potential_novel_variant"].sum()),
            "confidence_breakdown":   conf_summary,
            "top_hosts":              host_summary,
            "sequencing_technologies": tech_summary,
        })
    df = pd.DataFrame(rows)
    p = os.path.join(out_dir, "distances", "summary_per_subtype.csv")
    df.to_csv(p, index=False)
    log.info(f"  ✔ {p}")
    return df


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  MAIN                                                                    ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def main():
    args = parse_args()
    out = args.output_dir

    for sub in ("fasta", "aligned", "distances", "trees", "clusters", "plots"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)

    log = setup_logging(out, args.verbose)
    log.info("=" * 70)
    log.info("  Blastocystis SSU rRNA Pipeline v5 — broad-first acquisition")
    log.info(f"  Started: {datetime.now().isoformat()}")
    log.info(f"  Run mode: {'SKIP-DOWNLOAD' if args.skip_download else 'FULL ACQUISITION'}")
    if args.resume:
        log.info("  Mode: RESUME (skipping accessions already in manifest)")
    log.info(f"  Strategy: BROAD-FIRST (max {args.max_broad} records)")
    log.info(f"  ST validation: {'OFF' if args.no_st_validation else 'ON'}")
    log.info(f"  Narrow coverage check: {'ON' if args.narrow_coverage_check else 'OFF'}")
    log.info(f"  Cross-ST dedup: {'OFF' if args.no_cross_st_dedup else 'ON'}")
    log.info(f"  Keep unassigned: {'ON' if args.keep_unassigned else 'OFF'}")
    log.info("=" * 70)

    # ── Subtype selection ────────────────────────────────────────────────────
    if args.subtypes:
        subtypes = [s.strip().upper() for s in args.subtypes.split(",")]
    else:
        subtypes = list(SUBTYPES_DEFAULT)
        if not args.include_nmast:
            subtypes = [st for st in subtypes if not st.startswith("NMAST")]
        if not args.include_contested:
            subtypes = [st for st in subtypes
                        if SUBTYPE_REGISTRY.get(st, {}).get("status") != "contested"]
        if args.status_filter:
            allowed = {s.strip().lower() for s in args.status_filter.split(",")}
            subtypes = [st for st in subtypes
                        if SUBTYPE_REGISTRY.get(st, {}).get("status", "").lower() in allowed]

    log.info(f"\nSubtypes selected: {len(subtypes)}")
    for st in subtypes:
        meta = SUBTYPE_REGISTRY.get(st, {})
        log.info(f"  {st:<8} [{meta.get('status','unknown'):<12}] "
                 f"[{'human+' if meta.get('human') else 'animal'}]  "
                 f"{meta.get('note','')[:60]}")

    # ── Aligner detection ────────────────────────────────────────────────────
    aligner = args.aligner
    if aligner == "auto":
        aligner = detect_aligner()
    use_kmer = aligner == "kmer"
    log.info(f"\nAligner: {aligner}")
    if use_kmer:
        log.warning("  No external aligner — using k-mer Jaccard distance.\n"
                    "  Install mafft for publication-quality analysis.")

    # ╔════════════════════════════════════════════════════════════════════╗
    # ║  ACQUISITION (v5 broad-first)                                      ║
    # ╚════════════════════════════════════════════════════════════════════╝
    all_records = []
    all_manifest = []

    if not args.skip_download:
        init_entrez(args.email, args.api_key, log)

        # ── --resume: skip accessions already in the manifest ──────────
        already_fetched = set()
        if args.resume:
            existing_manifest_path = os.path.join(out, "fasta", "manifest.csv")
            if os.path.exists(existing_manifest_path):
                existing = pd.read_csv(existing_manifest_path)
                already_fetched = set(existing["accession"].tolist()) \
                                  | {a.split(".")[0] for a in existing["accession"]}
                log.info(f"\n[RESUME] {len(already_fetched)} accessions already in manifest;"
                         f" will skip them.")
                # Carry forward the existing manifest as a starting point
                for _, row in existing.iterrows():
                    all_manifest.append(row.to_dict())

        # ── STEP 1A — broad query (the primary acquisition route) ────────
        log.info("\n[STEP 1A] BROAD ACQUISITION — genus-wide Blastocystis SSU rRNA")
        broad_query = broad_blastocystis_query()
        broad_records_list = fetch_records_by_query(
            broad_query, args.max_broad, log,
            batch_size=args.fetch_batch_size,
            retries=args.fetch_retries,
        )
        # Filter out resume-already-have, and convert list → dict keyed by id
        broad_records = {r.id: r for r in broad_records_list
                          if r.id not in already_fetched
                          and r.id.split(".")[0] not in already_fetched}
        if args.resume and len(broad_records) < len(broad_records_list):
            log.info(f"  --resume: skipped {len(broad_records_list) - len(broad_records)}"
                     f" already-fetched records")
        log.info(f"  Broad query yielded {len(broad_records)} new records.")

        # ── STEP 1B — literature-curated reference accessions ────────────
        log.info("\n[STEP 1B] CURATED REFERENCE ACCESSIONS")
        ref_to_fetch = [a for a in REFERENCE_ST_MAP
                        if a not in already_fetched
                        and a not in broad_records
                        and a + ".1" not in broad_records]
        log.info(f"  {len(REFERENCE_ST_MAP)} known references, "
                 f"{len(ref_to_fetch)} not yet fetched.")
        ref_records = fetch_accessions(
            ref_to_fetch, log,
            batch_size=args.fetch_batch_size,
            retries=args.fetch_retries,
        ) if ref_to_fetch else {}

        # ── STEP 1C — optional narrow per-ST coverage check ──────────────
        narrow_records = {}
        if args.narrow_coverage_check:
            log.info("\n[STEP 1C] NARROW COVERAGE CHECK — per-ST title-scoped queries")
            for st in tqdm(subtypes, desc="  narrow check"):
                recs = fetch_records_by_query(
                    build_narrow_query(st),
                    args.max_per_st, log,
                    batch_size=args.fetch_batch_size,
                    retries=args.fetch_retries,
                )
                for r in recs:
                    if (r.id not in already_fetched
                            and r.id not in broad_records
                            and r.id not in ref_records
                            and r.id not in narrow_records):
                        narrow_records[r.id] = r
            log.info(f"  Narrow check found {len(narrow_records)} additional records "
                     f"not retrieved by the broad query.")
        else:
            log.info("\n[STEP 1C] Skipping narrow coverage check (use --narrow-coverage-check to enable)")

        # ── STEP 1D — record-quality and length filters ────────────────────
        # Quality filter: drop RefSeq predicted/curated records (these are
        # not deposited genomic SSU rRNA — they are computationally-predicted
        # transcripts/proteins that contaminate ST clusters with extreme
        # distance outliers when aligned against real SSU rRNA sequences).
        #
        # RefSeq accession prefixes we reject:
        #   XM_ : predicted mRNA model         XP_ : predicted protein model
        #   XR_ : predicted ncRNA model        YP_ : protein from genome project
        #   NM_ : curated mRNA                 NP_ : curated protein
        #   NR_ : curated ncRNA                WP_ : non-redundant protein
        # We also drop any record whose /mol_type is incompatible with a
        # deposited SSU rRNA gene (we accept "genomic DNA", "rRNA",
        # "transcribed RNA", "DNA"; we reject "mRNA" and "protein").
        REFSEQ_PREDICTED_PREFIXES = (
            "XM_", "XR_", "XP_",
            "NM_", "NR_", "NP_",
            "YP_", "WP_",
        )
        ACCEPTED_MOL_TYPES = {
            "", "genomic dna", "rrna", "trna", "dna",
            "transcribed rna", "other rna", "ss-rna", "ds-rna",
            "unassigned rna", "genomic rna",
        }

        def _is_refseq_predicted(rec) -> bool:
            return rec.id.startswith(REFSEQ_PREDICTED_PREFIXES)

        def _bad_mol_type(rec) -> str:
            """Return offending mol_type, or '' if record is acceptable."""
            for feat in rec.features:
                if feat.type == "source":
                    mt = (feat.qualifiers.get("mol_type", [""])[0] or "").strip().lower()
                    if mt and mt not in ACCEPTED_MOL_TYPES:
                        # Strong reject: clearly not an rRNA-gene deposit
                        if "protein" in mt or mt == "mrna":
                            return mt
                    return ""
            return ""

        def _quality_filter(d: dict, label: str) -> dict:
            kept = {}
            dropped_refseq = []
            dropped_moltype = []
            for k, rec in d.items():
                if _is_refseq_predicted(rec):
                    dropped_refseq.append(rec.id)
                    continue
                bad = _bad_mol_type(rec)
                if bad:
                    dropped_moltype.append((rec.id, bad))
                    continue
                kept[k] = rec
            if dropped_refseq:
                log.warning(f"  Quality filter on {label}: dropped "
                            f"{len(dropped_refseq)} RefSeq-predicted records "
                            f"(XM_/XR_/etc): {dropped_refseq[:5]}"
                            + ("..." if len(dropped_refseq) > 5 else ""))
            if dropped_moltype:
                examples = [f"{acc} (mol_type={mt})" for acc, mt in dropped_moltype[:3]]
                log.warning(f"  Quality filter on {label}: dropped "
                            f"{len(dropped_moltype)} records with non-rRNA "
                            f"mol_type: {examples}"
                            + ("..." if len(dropped_moltype) > 3 else ""))
            return kept

        def _len_filter(d: dict, label: str) -> dict:
            kept = {k: r for k, r in d.items()
                    if args.min_len <= len(r.seq) <= args.max_len}
            if len(kept) < len(d):
                log.info(f"  Length filter on {label}: {len(d)} → {len(kept)}")
            return kept

        log.info("\n[STEP 1D-i] Record-quality filtering (RefSeq-predicted, mol_type)")
        broad_records   = _quality_filter(broad_records,   "broad")
        ref_records     = _quality_filter(ref_records,     "reference")
        narrow_records  = _quality_filter(narrow_records,  "narrow check")

        log.info("\n[STEP 1D-ii] Length-filtering all batches")
        broad_records   = _len_filter(broad_records,   "broad")
        ref_records     = _len_filter(ref_records,     "reference")
        narrow_records  = _len_filter(narrow_records,  "narrow check")

        # ── STEP 1E — assign each record to an ST ────────────────────────
        log.info("\n[STEP 1E] ASSIGNING ST per record")
        log.info("  Two-stage: (i) parse source qualifiers; (ii) for unlabelled, "
                 "place against curated + high-confidence seeds.")

        retrieved_at = datetime.now(timezone.utc).isoformat()
        validate_st = not args.no_st_validation

        # Pool 1: all newly-fetched records, with their provenance
        # (we tag each with the acquisition source)
        all_new = []
        for rec_id, rec in broad_records.items():
            all_new.append((rec, "broad_query", broad_query))
        for rec_id, rec in ref_records.items():
            all_new.append((rec, "curated_reference", "REFERENCE_LOOKUP"))
        for rec_id, rec in narrow_records.items():
            all_new.append((rec, "narrow_coverage_check", "narrow per-ST query"))

        # ── 1E-i. Try to assign each record from its own source qualifiers ──
        # Records that declare a single ST → high confidence assignment
        # Records that declare nothing → mark for placement step
        labelled = []           # (rec, st, confidence, reason, source, query)
        needs_placement = []    # records lacking ST declaration

        for rec, source, query in all_new:
            # Curated reference: use the literature-mapped ST directly
            if source == "curated_reference":
                acc_bare = rec.id.split(".")[0]
                st = REFERENCE_ST_MAP.get(acc_bare) or REFERENCE_ST_MAP.get(rec.id)
                if st:
                    labelled.append((rec, st, "curated_reference",
                                     f"literature-curated reference (→ {st})",
                                     source, query))
                continue

            declared = extract_declared_sts(rec)
            if len(declared) == 1:
                st = next(iter(declared))
                if validate_st:
                    labelled.append((rec, st, "high",
                                     f"source declares {st} only", source, query))
                else:
                    labelled.append((rec, st, "unchecked",
                                     "validation disabled", source, query))
            elif len(declared) > 1:
                # Multiple STs declared — keep but flag as ambiguous
                # (any of them might be right). For broad-first we keep with the
                # most common ST, but flag medium.
                # If the record only declares two and one is mentioned in title, use that
                st = sorted(declared)[0]  # arbitrary but deterministic
                labelled.append((rec, st, "medium",
                                 f"source declares multiple STs {sorted(declared)}; "
                                 f"defaulting to {st}", source, query))
            else:
                needs_placement.append((rec, source, query))

        log.info(f"  {len(labelled)} records labelled by source qualifiers; "
                 f"{len(needs_placement)} need phylogenetic placement.")

        # ── 1E-ii. Build seed set from labelled records (high confidence only) ─
        seed_records = [rec for rec, st, conf, *_ in labelled
                        if conf in ("high", "curated_reference")]
        seed_st_map = {rec.id: st for rec, st, conf, *_ in labelled
                       if conf in ("high", "curated_reference")}

        # Also add any high-confidence records already in resume manifest
        if all_manifest:
            existing_high = {m["accession"]: m["subtype"] for m in all_manifest
                             if m.get("st_confidence") in ("high", "curated_reference")}
            # Resume records are already in manifest but we don't have SeqRecords
            # for them here (they're in the on-disk FASTAs). For placement we just
            # need their sequences. Read from per-ST FASTAs if available.
            for st_fasta in os.listdir(os.path.join(out, "fasta")):
                if not st_fasta.endswith("_filtered.fasta"):
                    continue
                p = os.path.join(out, "fasta", st_fasta)
                for rec in SeqIO.parse(p, "fasta"):
                    if rec.id in existing_high and rec.id not in seed_st_map:
                        seed_records.append(rec)
                        seed_st_map[rec.id] = existing_high[rec.id]

        log.info(f"  Placement seed set: {len(seed_records)} sequences "
                 f"spanning {len(set(seed_st_map.values()))} STs.")

        # ── 1E-iii. Place unlabelled records ──────────────────────────────
        placed = []           # (rec, st, confidence, reason, source, query)
        unassigned = []       # (rec, source, query) — could not be placed

        if needs_placement and seed_records:
            placement_results = place_against_seeds(
                [r for r, _, _ in needs_placement],
                seed_records, seed_st_map, log,
                max_dist=args.max_placement_dist,
            )
            for rec, source, query in needs_placement:
                st, dist, conf, reason = placement_results.get(
                    rec.id, ("Unplaced", float("inf"), "rejected", "no placement"))
                if st == "Unplaced":
                    unassigned.append((rec, source, query))
                else:
                    placed.append((rec, st, conf, reason, source, query))
        else:
            # No seeds (first run with no existing data) — everything is unassigned
            unassigned = needs_placement

        log.info(f"  Placement: {len(placed)} placed, {len(unassigned)} unassigned.")
        if placed:
            placement_confs = Counter(c for _, _, c, _, _, _ in placed)
            log.info("    Placement confidence breakdown:")
            for c, n in placement_confs.most_common():
                log.info(f"      {c:<10}: {n}")

        # ── 1E-iv. Build the new records + manifest rows ──────────────────
        new_records = []
        new_manifest = []

        def _build_row(rec, st, confidence, reason, source, query):
            meta = extract_metadata_block(rec)
            fasta_rec = SeqRecord(
                rec.seq, id=rec.id,
                description=(
                    f"{rec.description} [subtype={st}] "
                    f"[host={meta['host']}] [country={meta['country']}] "
                    f"[date={meta['collection_date']}] "
                    f"[isolation_source={meta['isolation_source']}] "
                    f"[seq_tech={meta['seq_tech']}] "
                    f"[confidence={confidence}] [host_source={meta['host_source']}]"
                )
            )
            row = {
                "accession":         rec.id,
                "subtype":           st,
                "length":            meta["length"],
                "md5":               meta["md5"],
                "host":              meta["host"],
                "host_source":       meta["host_source"],
                "country":           meta["country"],
                "collection_date":   meta["collection_date"],
                "isolation_source":  meta["isolation_source"],
                "seq_tech":          meta["seq_tech"],
                "st_confidence":     confidence,
                "st_validation":     reason,
                "acquisition_pass":  source,
                "retrieved_at_utc":  retrieved_at,
                "ncbi_query":        query,
            }
            return fasta_rec, row

        for rec, st, conf, reason, source, query in labelled + placed:
            fasta_rec, row = _build_row(rec, st, conf, reason, source, query)
            new_records.append(fasta_rec)
            new_manifest.append(row)

        # ── 1E-iv-bis. Novel-subtype clustering of unassigned records ────
        # Any record that placement could not assign to a known ST is a
        # candidate for novel-subtype detection. We cluster these among
        # themselves and label coherent groups novel_1, novel_2, ...
        novel_assignments: dict = {}
        if args.cluster_novel_subtypes and unassigned:
            log.info("\n  Novel-subtype clustering of unplaced records")
            novel_assignments = cluster_novel_subtypes(
                [r for r, _, _ in unassigned], log,
                cluster_distance_cutoff=args.novel_cluster_cutoff,
                min_cluster_size=args.novel_min_size,
            )

        if args.keep_unassigned:
            for rec, source, query in unassigned:
                novel_label = novel_assignments.get(rec.id)
                if novel_label is not None:
                    fasta_rec, row = _build_row(
                        rec, novel_label, "novel_cluster",
                        f"clustered with {sum(1 for v in novel_assignments.values() if v == novel_label)} "
                        f"other unplaced records at "
                        f"Jaccard d<{args.novel_cluster_cutoff}",
                        source, query,
                    )
                else:
                    fasta_rec, row = _build_row(
                        rec, "ST_unassigned", "unplaced",
                        "no source ST tag, >max-placement-dist from all seeds, "
                        "and no coherent novel cluster",
                        source, query,
                    )
                new_records.append(fasta_rec)
                new_manifest.append(row)

            n_in_novel = len(novel_assignments)
            n_singletons = len(unassigned) - n_in_novel
            log.info(f"  Of {len(unassigned)} unplaced records: "
                     f"{n_in_novel} assigned to provisional novel clusters, "
                     f"{n_singletons} remain ST_unassigned.")

        log.info(f"\n  Acquired {len(new_records)} new records this run "
                 f"(labelled: {len(labelled)}, placed: {len(placed)}, "
                 f"novel: {len(novel_assignments)}, "
                 f"unassigned: {len(unassigned) - len(novel_assignments)}).")

        # ── STEP 1F — Merge with resume manifest, cross-ST dedup ──────────
        # If resuming, we have all_manifest from disk; need to also load the
        # SeqRecords for the resumed accessions to write the combined output.
        if args.resume and all_manifest:
            log.info("\n[STEP 1F] Loading resumed records from disk")
            existing_recs = []
            for st_fasta in os.listdir(os.path.join(out, "fasta")):
                if st_fasta.endswith("_filtered.fasta"):
                    existing_recs.extend(SeqIO.parse(
                        os.path.join(out, "fasta", st_fasta), "fasta"))
            existing_ids = {r.id for r in existing_recs}
            # Drop any that are also in new_records (newest wins)
            new_ids = {r.id for r in new_records}
            existing_recs = [r for r in existing_recs if r.id not in new_ids]
            all_records = existing_recs + new_records
            log.info(f"  {len(existing_recs)} from previous run + "
                     f"{len(new_records)} new = {len(all_records)} total.")
        else:
            all_records = new_records
            all_manifest = new_manifest

        # If resuming, all_manifest already has rows; we need to append new ones
        if args.resume:
            existing_accs = {m["accession"] for m in all_manifest}
            for row in new_manifest:
                if row["accession"] not in existing_accs:
                    all_manifest.append(row)

        # ── STEP 1G — Cross-ST dedup ─────────────────────────────────────
        if not args.no_cross_st_dedup and len(all_records) > 1:
            log.info("\n[STEP 1G] Cross-ST deduplication")
            kept = deduplicate_across_sts(all_records, all_manifest, log)
            kept_ids = {r.id for r in kept}
            all_records = kept
            all_manifest = [m for m in all_manifest if m["accession"] in kept_ids]

        # ── STEP 1H — Write per-ST FASTAs + manifest ─────────────────────
        log.info("\n[STEP 1H] Writing FASTAs and manifest")
        # Clear old per-ST FASTAs to avoid stale carry-over
        for f in os.listdir(os.path.join(out, "fasta")):
            if f.endswith("_filtered.fasta"):
                os.remove(os.path.join(out, "fasta", f))

        by_st = defaultdict(list)
        st_of = {m["accession"]: m["subtype"] for m in all_manifest}
        for r in all_records:
            by_st[st_of.get(r.id, "ST_unassigned")].append(r)
        for st, recs in by_st.items():
            p = os.path.join(out, "fasta", f"{st}_filtered.fasta")
            SeqIO.write(recs, p, "fasta")
            log.info(f"  {st}: {len(recs)} → {p}")
        combined = os.path.join(out, "fasta", "all_subtypes.fasta")
        SeqIO.write(all_records, combined, "fasta")
        manifest_df = pd.DataFrame(all_manifest)
        manifest_df.to_csv(os.path.join(out, "fasta", "manifest.csv"), index=False)
        log.info(f"  Combined: {len(all_records)} → {combined}")
        log.info(f"  Manifest: {os.path.join(out, 'fasta', 'manifest.csv')}")

    else:
        log.info("\n[STEP 1] Skipping download — loading existing FASTAs")
        combined = os.path.join(out, "fasta", "all_subtypes.fasta")
        all_records = list(SeqIO.parse(combined, "fasta"))
        log.info(f"  Loaded {len(all_records)} sequences from {combined}")
        manifest_path = os.path.join(out, "fasta", "manifest.csv")
        if os.path.exists(manifest_path):
            manifest_df = pd.read_csv(manifest_path)
        else:
            manifest_df = pd.DataFrame()

    if not all_records:
        log.error("No sequences available. Exiting.")
        sys.exit(1)

    # ╔════════════════════════════════════════════════════════════════════╗
    # ║  STEP 2 — ALIGNMENT                                                ║
    # ╚════════════════════════════════════════════════════════════════════╝
    log.info(f"\n[STEP 2] Aligning ({aligner})...")
    aligned_path = os.path.join(out, "aligned", "all_subtypes_aligned.fasta")
    align_sequences(combined, aligned_path, aligner, args.threads, log)
    aligned_records = list(SeqIO.parse(aligned_path, "fasta")) or all_records
    log.info(f"  {len(aligned_records)} sequences ready for analysis")

    # ╔════════════════════════════════════════════════════════════════════╗
    # ║  STEP 3 — DISTANCES                                                ║
    # ╚════════════════════════════════════════════════════════════════════╝
    log.info("\n[STEP 3] Pairwise distance matrix")
    D, ids = compute_distance_matrix(aligned_records, use_kmer, args.jc_correct, log)
    pd.DataFrame(D, index=ids, columns=ids).to_csv(
        os.path.join(out, "distances", "distance_matrix.csv"), float_format="%.5f")
    np.savez_compressed(
        os.path.join(out, "distances", "distance_matrix.npz"),
        D=D, ids=np.array(ids))
    log.info(f"  ✔ distance_matrix.csv + .npz")

    # ╔════════════════════════════════════════════════════════════════════╗
    # ║  STEP 4 — CLUSTERING                                               ║
    # ╚════════════════════════════════════════════════════════════════════╝
    log.info("\n[STEP 4] Clustering")
    cluster_df, Z = cluster_sequences(
        D, ids, aligned_records, log,
        thresh_cluster=args.thresh_cluster,
        thresh_subtype=args.thresh_subtype,
        thresh_variant=args.thresh_variant,
        novel_distance_multiplier=args.novel_distance_multiplier,
    )
    cluster_path = os.path.join(out, "clusters", "cluster_assignments.csv")
    cluster_df.to_csv(cluster_path, index=False)
    log.info(f"  ✔ {cluster_path}")

    # Save run config
    with open(os.path.join(out, "run_config.json"), "w") as f:
        json.dump({
            "version": "v5",
            "args": vars(args),
            "aligner": aligner,
            "n_sequences": int(len(aligned_records)),
            "n_subtypes": int(cluster_df["subtype"].nunique()),
            "started": datetime.now().isoformat(),
        }, f, indent=2, default=str)

    # ╔════════════════════════════════════════════════════════════════════╗
    # ║  STEP 5 — NJ TREE                                                  ║
    # ╚════════════════════════════════════════════════════════════════════╝
    log.info("\n[STEP 5] Neighbor-joining tree")
    if use_kmer:
        log.warning("  Skipping NJ: k-mer Jaccard not additive — invalid for NJ.")
    elif len(aligned_records) <= 1000:
        try:
            newick = neighbor_joining(D, ids)
            tree_p = os.path.join(out, "trees", "nj_tree.nwk")
            with open(tree_p, "w") as f:
                f.write(newick)
            log.info(f"  ✔ {tree_p}")
        except Exception as e:
            log.warning(f"  NJ failed: {e}")
    else:
        log.info(f"  Skipping NJ (>1000 seqs). For publication-grade tree:")
        log.info(f"    iqtree -s {aligned_path} -m GTR+G -bb 1000 -nt AUTO")

    # ╔════════════════════════════════════════════════════════════════════╗
    # ║  STEP 6 — ORDINATION                                               ║
    # ╚════════════════════════════════════════════════════════════════════╝
    log.info("\n[STEP 6] Ordination")
    coords, method, var = compute_embedding(D, log, seed=args.seed)

    # ╔════════════════════════════════════════════════════════════════════╗
    # ║  STEP 7 — PLOTS & REPORTS                                          ║
    # ╚════════════════════════════════════════════════════════════════════╝
    log.info("\n[STEP 7] Plots & reports")
    plot_ordination(coords, cluster_df, method, var, out, log)
    plot_distance_heatmap(D, cluster_df, out, log)
    plot_intra_inter(D, cluster_df, out, log, thresh_subtype=args.thresh_subtype)
    plot_technology_breakdown(cluster_df, out, log)
    plot_host_distribution(cluster_df, out, log)
    plot_geographic_temporal(cluster_df, out, log)
    plot_novel_clusters(cluster_df, out, log)
    if not manifest_df.empty:
        plot_acquisition_breakdown(manifest_df, out, log)
    summary_df = write_summary_table(D, cluster_df, out, log)
    build_interactive_report(coords, cluster_df, D, method, out, log)

    # ╔════════════════════════════════════════════════════════════════════╗
    # ║  FINAL SUMMARY                                                     ║
    # ╚════════════════════════════════════════════════════════════════════╝
    log.info("\n" + "=" * 70)
    log.info("  PIPELINE COMPLETE — v5")
    log.info("=" * 70)
    log.info(f"  Sequences analysed       : {len(aligned_records)}")
    log.info(f"  Subtypes in query        : {len(subtypes)}")
    log.info(f"  Subtypes with sequences  : {cluster_df['subtype'].nunique()}")
    log.info(f"  Variant-level clusters   : {cluster_df['cluster_variant'].nunique()}")
    log.info(f"  Novel-variant candidates : {int(cluster_df['potential_novel_variant'].sum())}")

    # Novel-subtype provisional clusters (v5 contribution)
    novel_sts = sorted(s for s in cluster_df["subtype"].unique()
                       if isinstance(s, str) and s.startswith("novel_"))
    n_unassigned = (cluster_df["subtype"] == "ST_unassigned").sum() \
                   if "subtype" in cluster_df.columns else 0
    if novel_sts or n_unassigned:
        log.info("\n  Provisional novel-subtype clusters:")
        for st in novel_sts:
            n = (cluster_df["subtype"] == st).sum()
            sample_hosts = (cluster_df[cluster_df["subtype"] == st]["host"]
                            .value_counts().head(2).index.tolist())
            log.info(f"    {st:<14}: n={n:<4}  top host(s): {sample_hosts}")
        if n_unassigned:
            log.info(f"    ST_unassigned : n={n_unassigned}  "
                     f"(singletons that did not form a coherent novel cluster)")

    if "st_confidence" in cluster_df.columns:
        log.info("\n  ST-claim confidence:")
        for c, n in cluster_df["st_confidence"].value_counts().items():
            log.info(f"    {c:<22}: {n}  ({100*n/len(cluster_df):.1f}%)")

    if "host_source" in cluster_df.columns:
        log.info("\n  Host derivation:")
        for s, n in cluster_df["host_source"].value_counts().items():
            log.info(f"    {s:<22}: {n}  ({100*n/len(cluster_df):.1f}%)")

    if "seq_technology" in cluster_df.columns:
        log.info("\n  Sequencing technologies:")
        for t, n in cluster_df["seq_technology"].value_counts().items():
            log.info(f"    {t:<22}: {n}  ({100*n/len(cluster_df):.1f}%)")

    log.info("\n  Outputs:")
    log.info(f"    Aligned FASTA    : {aligned_path}")
    log.info(f"    Distance matrix  : distances/distance_matrix.{{csv,npz}}")
    log.info(f"    Cluster table    : {cluster_path}")
    log.info(f"    Per-ST summary   : distances/summary_per_subtype.csv")
    log.info(f"    Manifest         : fasta/manifest.csv")
    log.info(f"    Run config       : run_config.json")
    log.info(f"    Plots            : plots/*.png")
    log.info(f"    Interactive      : plots/interactive_report.html")
    log.info("\n  For publication-quality phylogeny:")
    log.info(f"    iqtree -s {aligned_path} -m GTR+G -bb 1000 -nt AUTO")
    log.info("=" * 70)

    if not summary_df.empty:
        log.info("\n  Per-subtype summary:")
        log.info(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
