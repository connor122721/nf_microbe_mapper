#!/usr/bin/env python3
"""
make_sample_report.py
---------------------
Render a single self-contained HTML report for one sample by parsing the
flagstat / idxstats / kraken2 / bracken outputs that Nextflow staged into
``--inputs-dir``.

Files are auto-detected by suffix:

    *.flagstat.txt            -> samtools flagstat (per reference)
    *.idxstats.txt            -> samtools idxstats (per reference)
    *.kraken2.report.txt      -> kraken2 report
    *.bracken.tsv             -> bracken abundance estimates

Stdlib only - no jinja, pandas, etc.  Container is python:3.10.
"""

from __future__ import annotations

import argparse
import csv
import html
import re
from datetime import datetime
from pathlib import Path
from typing import Iterable


# -----------------------------------------------------------------------------
# Parsers
# -----------------------------------------------------------------------------

FLAGSTAT_KEYS = (
    ("total",         re.compile(r"^(\d+)\s*\+\s*\d+\s+in total")),
    ("mapped",        re.compile(r"^(\d+)\s*\+\s*\d+\s+mapped")),
    ("primary",       re.compile(r"^(\d+)\s*\+\s*\d+\s+primary$")),
    ("secondary",     re.compile(r"^(\d+)\s*\+\s*\d+\s+secondary")),
    ("supplementary", re.compile(r"^(\d+)\s*\+\s*\d+\s+supplementary")),
    ("duplicates",    re.compile(r"^(\d+)\s*\+\s*\d+\s+duplicates")),
)


def parse_flagstat(path: Path) -> dict[str, int]:
    out: dict[str, int] = {k: 0 for k, _ in FLAGSTAT_KEYS}
    for line in path.read_text().splitlines():
        for key, pat in FLAGSTAT_KEYS:
            m = pat.match(line.strip())
            if m:
                out[key] = int(m.group(1))
    return out


def parse_idxstats(path: Path) -> list[tuple[str, int, int, int]]:
    """Returns list of (chrom, length, mapped, unmapped)."""
    rows: list[tuple[str, int, int, int]] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        try:
            rows.append((parts[0], int(parts[1]), int(parts[2]), int(parts[3])))
        except ValueError:
            continue
    return rows


def parse_kraken_report(path: Path, top_n: int = 15) -> list[dict]:
    """
    Kraken2 report columns:
      0=pct, 1=reads_clade, 2=reads_taxon, 3=rank, 4=taxid, 5=name
      (extra columns when --report-minimizer-data is on; we ignore them)
    Returns the top_n entries at species rank, then falls back to genus, etc.
    """
    rows: list[dict] = []
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                pct = float(parts[0])
            except ValueError:
                continue
            rows.append({
                "pct":          pct,
                "reads_clade":  int(parts[1]),
                "reads_taxon":  int(parts[2]),
                "rank":         parts[3].strip(),
                "taxid":        parts[4].strip(),
                "name":         parts[5].strip(),
            })

    # Prefer species; if there are very few, blend in genus
    species = [r for r in rows if r["rank"] == "S"]
    species.sort(key=lambda r: r["pct"], reverse=True)
    if len(species) >= 5:
        return species[:top_n]
    genus = [r for r in rows if r["rank"] == "G"]
    genus.sort(key=lambda r: r["pct"], reverse=True)
    return (species + genus)[:top_n]


def parse_bracken(path: Path, top_n: int = 15) -> list[dict]:
    """
    Bracken TSV columns:
      name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads,
      added_reads, new_est_reads, fraction_total_reads
    """
    rows: list[dict] = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                rows.append({
                    "name":           row["name"],
                    "taxid":          row.get("taxonomy_id", ""),
                    "level":          row.get("taxonomy_lvl", ""),
                    "kraken_reads":   int(row.get("kraken_assigned_reads", 0) or 0),
                    "new_est_reads":  int(row.get("new_est_reads", 0) or 0),
                    "fraction":       float(row.get("fraction_total_reads", 0) or 0),
                })
            except (KeyError, ValueError):
                continue
    rows.sort(key=lambda r: r["fraction"], reverse=True)
    return rows[:top_n]


# -----------------------------------------------------------------------------
# HTML rendering (no external deps)
# -----------------------------------------------------------------------------

CSS = """
body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif;
     margin:0;padding:24px;color:#1f2328;background:#fafbfc;max-width:1100px;margin:auto;}
h1{margin-top:0;border-bottom:2px solid #d0d7de;padding-bottom:8px;}
h2{margin-top:32px;color:#0969da;}
.meta{color:#57606a;font-size:13px;margin-bottom:16px;}
table{border-collapse:collapse;width:100%;margin-top:8px;background:#fff;
      box-shadow:0 1px 0 rgba(27,31,36,0.04);}
th,td{padding:8px 10px;border-bottom:1px solid #d0d7de;text-align:left;font-size:13px;}
th{background:#f6f8fa;font-weight:600;}
tr:hover td{background:#f6f8fa;}
.num{text-align:right;font-variant-numeric:tabular-nums;}
.bar{background:#0969da;height:10px;border-radius:5px;}
.bar-bg{background:#eaeef2;border-radius:5px;width:160px;display:inline-block;vertical-align:middle;}
.tag{display:inline-block;padding:2px 8px;border-radius:12px;font-size:11px;
     background:#ddf4ff;color:#0969da;margin-right:6px;}
.callout{background:#fff8c5;border:1px solid #d4a72c;padding:8px 12px;border-radius:6px;
         font-size:13px;}
footer{margin-top:48px;color:#57606a;font-size:12px;text-align:center;}
"""


def fmt_int(n: int) -> str:
    return f"{n:,}"


def pct_bar(value: float, max_val: float = 100.0) -> str:
    width = max(1, min(160, int(160 * value / max_val))) if max_val else 1
    return f'<span class="bar-bg"><span class="bar" style="width:{width}px;display:block"></span></span>'


def render_html(sample: str,
                flagstats: list[tuple[str, dict[str, int]]],
                idxstats: list[tuple[str, list]],
                kraken_top: list[dict],
                bracken_top: list[dict]) -> str:

    parts: list[str] = []
    parts.append(f"<!doctype html><html><head><meta charset='utf-8'>"
                 f"<title>{html.escape(sample)} - Sample Report</title>"
                 f"<style>{CSS}</style></head><body>")
    parts.append(f"<h1>Sample report: <code>{html.escape(sample)}</code></h1>")
    parts.append(f"<div class='meta'>Generated {datetime.utcnow().isoformat(timespec='seconds')}Z "
                 f"by <code>pb-hifi-mapping-nf</code></div>")

    # ---------- mapping summary ----------
    parts.append("<h2>Mapping summary (per reference)</h2>")
    if not flagstats:
        parts.append("<div class='callout'>No flagstat files were produced.</div>")
    else:
        parts.append("<table><thead><tr>"
                     "<th>Reference</th>"
                     "<th class='num'>Total reads</th>"
                     "<th class='num'>Mapped</th>"
                     "<th class='num'>% Mapped</th>"
                     "<th class='num'>Primary</th>"
                     "<th class='num'>Secondary</th>"
                     "<th class='num'>Supplementary</th>"
                     "<th class='num'>Duplicates</th>"
                     "</tr></thead><tbody>")
        for ref_id, fs in flagstats:
            total  = fs["total"] or 0
            mapped = fs["mapped"] or 0
            pct    = (mapped / total * 100.0) if total else 0.0
            parts.append(
                f"<tr>"
                f"<td><code>{html.escape(ref_id)}</code></td>"
                f"<td class='num'>{fmt_int(total)}</td>"
                f"<td class='num'>{fmt_int(mapped)}</td>"
                f"<td class='num'>{pct:.2f}% {pct_bar(pct)}</td>"
                f"<td class='num'>{fmt_int(fs['primary'])}</td>"
                f"<td class='num'>{fmt_int(fs['secondary'])}</td>"
                f"<td class='num'>{fmt_int(fs['supplementary'])}</td>"
                f"<td class='num'>{fmt_int(fs['duplicates'])}</td>"
                f"</tr>"
            )
        parts.append("</tbody></table>")

    # ---------- per-contig top hits ----------
    parts.append("<h2>Top contigs by mapped reads</h2>")
    any_idx = False
    for ref_id, rows in idxstats:
        rows_sorted = sorted([r for r in rows if r[0] != "*"], key=lambda r: r[2], reverse=True)[:10]
        if not rows_sorted:
            continue
        any_idx = True
        parts.append(f"<h3 style='margin-top:24px;font-size:15px'>Reference: "
                     f"<code>{html.escape(ref_id)}</code></h3>")
        parts.append("<table><thead><tr>"
                     "<th>Contig</th>"
                     "<th class='num'>Length (bp)</th>"
                     "<th class='num'>Mapped reads</th>"
                     "<th class='num'>Unmapped (this contig)</th>"
                     "</tr></thead><tbody>")
        for chrom, length, mapped, unmapped in rows_sorted:
            parts.append(
                f"<tr>"
                f"<td><code>{html.escape(chrom)}</code></td>"
                f"<td class='num'>{fmt_int(length)}</td>"
                f"<td class='num'>{fmt_int(mapped)}</td>"
                f"<td class='num'>{fmt_int(unmapped)}</td>"
                f"</tr>"
            )
        parts.append("</tbody></table>")
    if not any_idx:
        parts.append("<div class='callout'>No idxstats data found.</div>")

    # ---------- taxonomy: Bracken first if present, else Kraken2 ----------
    parts.append("<h2>Taxonomic profile</h2>")
    if bracken_top:
        parts.append(f"<div class='meta'>Source: <span class='tag'>Bracken</span> "
                     f"(top {len(bracken_top)} taxa, abundance-corrected)</div>")
        parts.append("<table><thead><tr>"
                     "<th>Taxon</th>"
                     "<th class='num'>Level</th>"
                     "<th class='num'>Kraken reads</th>"
                     "<th class='num'>Bracken est. reads</th>"
                     "<th class='num'>% of sample</th>"
                     "</tr></thead><tbody>")
        max_frac = max((b["fraction"] for b in bracken_top), default=1.0) or 1.0
        for b in bracken_top:
            parts.append(
                f"<tr>"
                f"<td>{html.escape(b['name'])}</td>"
                f"<td class='num'>{html.escape(b['level'])}</td>"
                f"<td class='num'>{fmt_int(b['kraken_reads'])}</td>"
                f"<td class='num'>{fmt_int(b['new_est_reads'])}</td>"
                f"<td class='num'>{b['fraction']*100:.3f}% {pct_bar(b['fraction']*100, max_frac*100)}</td>"
                f"</tr>"
            )
        parts.append("</tbody></table>")
    elif kraken_top:
        parts.append(f"<div class='meta'>Source: <span class='tag'>Kraken2</span> "
                     f"(top {len(kraken_top)} taxa)</div>")
        parts.append("<table><thead><tr>"
                     "<th>Taxon</th>"
                     "<th class='num'>Rank</th>"
                     "<th class='num'>Reads (clade)</th>"
                     "<th class='num'>% of sample</th>"
                     "</tr></thead><tbody>")
        max_pct = max((k["pct"] for k in kraken_top), default=1.0) or 1.0
        for k in kraken_top:
            parts.append(
                f"<tr>"
                f"<td>{html.escape(k['name'])}</td>"
                f"<td class='num'>{html.escape(k['rank'])}</td>"
                f"<td class='num'>{fmt_int(k['reads_clade'])}</td>"
                f"<td class='num'>{k['pct']:.3f}% {pct_bar(k['pct'], max_pct)}</td>"
                f"</tr>"
            )
        parts.append("</tbody></table>")
    else:
        parts.append("<div class='callout'>Taxonomic profiling was not run for this sample "
                     "(no <code>--kraken2_db</code> provided).</div>")

    parts.append("<footer>pb-hifi-mapping-nf &middot; Connor S. Murray, PhD</footer>")
    parts.append("</body></html>")
    return "\n".join(parts)


# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------

def discover(inputs_dir: Path) -> dict[str, list[Path]]:
    bucket: dict[str, list[Path]] = {
        "flagstat": [], "idxstats": [], "kraken_report": [], "bracken": [],
    }
    for f in sorted(inputs_dir.iterdir()):
        if not f.is_file():
            continue
        name = f.name.lower()
        if name.endswith(".flagstat.txt"):
            bucket["flagstat"].append(f)
        elif name.endswith(".idxstats.txt"):
            bucket["idxstats"].append(f)
        elif name.endswith(".kraken2.report.txt"):
            bucket["kraken_report"].append(f)
        elif name.endswith(".bracken.tsv"):
            bucket["bracken"].append(f)
    return bucket


def ref_id_from(path: Path, suffix: str) -> str:
    # foo.vs.bar.flagstat.txt  -> bar
    n = path.name
    if n.lower().endswith(suffix):
        n = n[: -len(suffix)]
    if ".vs." in n:
        return n.split(".vs.", 1)[1]
    return n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample",      required=True)
    ap.add_argument("--inputs-dir",  required=True, type=Path)
    ap.add_argument("--html-out",    required=True, type=Path)
    ap.add_argument("--tsv-out",     required=True, type=Path)
    args = ap.parse_args()

    files = discover(args.inputs_dir)

    flagstats = [(ref_id_from(p, ".flagstat.txt"), parse_flagstat(p)) for p in files["flagstat"]]
    idxstats  = [(ref_id_from(p, ".idxstats.txt"), parse_idxstats(p))  for p in files["idxstats"]]
    kraken_top  = parse_kraken_report(files["kraken_report"][0]) if files["kraken_report"] else []
    bracken_top = parse_bracken(files["bracken"][0])             if files["bracken"]       else []

    # HTML
    args.html_out.write_text(render_html(args.sample, flagstats, idxstats, kraken_top, bracken_top))

    # Flat TSV summary - one row per reference
    with args.tsv_out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "reference", "total_reads", "mapped_reads", "pct_mapped",
                    "primary", "secondary", "supplementary", "duplicates",
                    "top_taxon", "top_taxon_pct"])
        top_name = ""
        top_pct  = ""
        if bracken_top:
            top_name = bracken_top[0]["name"]
            top_pct  = f"{bracken_top[0]['fraction']*100:.3f}"
        elif kraken_top:
            top_name = kraken_top[0]["name"]
            top_pct  = f"{kraken_top[0]['pct']:.3f}"
        for ref_id, fs in flagstats:
            total  = fs["total"] or 0
            mapped = fs["mapped"] or 0
            pct    = (mapped / total * 100.0) if total else 0.0
            w.writerow([args.sample, ref_id, total, mapped, f"{pct:.4f}",
                        fs["primary"], fs["secondary"], fs["supplementary"], fs["duplicates"],
                        top_name, top_pct])
        if not flagstats:
            w.writerow([args.sample, "", 0, 0, "0.0000", 0, 0, 0, 0, top_name, top_pct])

    print(f"OK: wrote {args.html_out} and {args.tsv_out}")


if __name__ == "__main__":
    main()
