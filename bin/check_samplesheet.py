#!/usr/bin/env python3
"""
check_samplesheet.py
--------------------
Standalone validator you can run BEFORE launching the pipeline:

    python bin/check_samplesheet.py samplesheet.csv

Verifies:
  - required columns present (sample, fastq)
  - sample IDs are unique
  - each fastq path resolves and has a recognized extension
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

REQUIRED_COLS = ("sample", "fastq")
OK_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".bam")


def fail(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samplesheet", help="Path to samplesheet CSV")
    ap.add_argument("--strict-paths", action="store_true",
                    help="Fail if any fastq path does not resolve on disk (default: warn)")
    args = ap.parse_args()

    sheet = Path(args.samplesheet)
    if not sheet.is_file():
        fail(f"samplesheet not found: {sheet}")

    with sheet.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            fail("samplesheet appears to be empty")
        missing = [c for c in REQUIRED_COLS if c not in reader.fieldnames]
        if missing:
            fail(f"samplesheet missing required column(s): {missing}.  Got: {reader.fieldnames}")

        seen: set[str] = set()
        problems: list[str] = []
        rows = 0

        for i, row in enumerate(reader, start=2):  # start=2 because header is line 1
            rows += 1
            sample = (row.get("sample") or "").strip()
            fastq = (row.get("fastq") or "").strip()

            if not sample:
                problems.append(f"line {i}: empty 'sample'")
            if not fastq:
                problems.append(f"line {i}: empty 'fastq'")
            if sample in seen:
                problems.append(f"line {i}: duplicate sample '{sample}'")
            seen.add(sample)

            if fastq and not any(fastq.lower().endswith(e) for e in OK_EXTS):
                problems.append(f"line {i}: '{fastq}' has unrecognized extension (allowed: {OK_EXTS})")

            if fastq:
                p = Path(fastq).expanduser()
                if not p.exists():
                    msg = f"line {i}: file does not exist: {fastq}"
                    if args.strict_paths:
                        problems.append(msg)
                    else:
                        print(f"WARN: {msg}", file=sys.stderr)

    if rows == 0:
        fail("samplesheet has no data rows")
    if problems:
        for p in problems:
            print(f"ERROR: {p}", file=sys.stderr)
        sys.exit(1)

    print(f"OK: {rows} sample(s) validated.")


if __name__ == "__main__":
    main()
