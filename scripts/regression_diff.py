#!/usr/bin/env python3
"""Diff two regression-summary CSVs and emit a markdown report.

Usage:
    regression_diff.py BASELINE_DIR NEW_DIR [--out OUT.md]

Each directory contains one or more `*.csv` files dumped by
`TransmissionLinelist.save_regression_summary`. CSVs sharing a basename
across the two directories are paired; CSVs present in one but not the
other are flagged. For each paired CSV, the script walks the
(section, quantity) rows and applies tolerances:

  - headline/rt:  |Δmedian| > max(0.02 * |baseline|, 0.1 * (hi - lo))
                  flags ❌; within ⅔ of that bound ⚠️; else ✅.
  - diagnostics:  rhat_max > 1.05 or grew ≥ 0.01 → ❌; ess_min dropped
                  > 30% or fell below 200 → ⚠️; divergences > 5 or
                  grew by ≥ 5 → ⚠️/❌ depending on absolute level.

The exit status is always 0 — the report is informational so it
surfaces in the PR comment without blocking the build. Treat the ❌
rows as the actionable signal.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path

REL_TOL = 0.02
CI_FRACTION = 0.1
WARN_FRACTION = 2.0 / 3.0


@dataclass(frozen=True)
class Row:
    section: str
    quantity: str
    median: float
    lower_95: float
    upper_95: float


def _read(path: Path) -> dict[tuple[str, str], Row]:
    out: dict[tuple[str, str], Row] = {}
    with path.open(newline="") as fh:
        for r in csv.DictReader(fh):
            section = r["section"]
            quantity = r["quantity"]
            row = Row(
                section,
                quantity,
                float(r["median"]),
                float(r["lower_95"]),
                float(r["upper_95"]),
            )
            out[(section, quantity)] = row
    return out


def _fmt(x: float, sig: int = 4) -> str:
    if math.isnan(x):
        return "—"
    if x == 0:
        return "0"
    return f"{x:#.{sig}g}"


def _interval(r: Row) -> str:
    if math.isnan(r.lower_95) or math.isnan(r.upper_95):
        return _fmt(r.median)
    return f"{_fmt(r.median)} [{_fmt(r.lower_95)}, {_fmt(r.upper_95)}]"


def _status_headline(base: Row, new: Row) -> tuple[str, str]:
    width = base.upper_95 - base.lower_95
    tol = max(REL_TOL * abs(base.median), CI_FRACTION * width)
    diff = abs(new.median - base.median)
    warn = WARN_FRACTION * tol
    if diff > tol:
        return ("FAIL", f"|Δ| {_fmt(diff)} > tol {_fmt(tol)}")
    if diff > warn:
        return ("WARN", f"|Δ| {_fmt(diff)} ≈ tol {_fmt(tol)}")
    return ("PASS", f"|Δ| {_fmt(diff)} ≤ {_fmt(warn)}")


def _status_diagnostic(base: Row, new: Row) -> tuple[str, str]:
    q = base.quantity
    if q == "rhat_max":
        if new.median > 1.05 or new.median - base.median >= 0.01:
            return ("FAIL", f"rhat_max {_fmt(new.median)} (was {_fmt(base.median)})")
        return ("PASS", f"rhat_max {_fmt(new.median)}")
    if q == "ess_min":
        if new.median < 200 or new.median < 0.7 * base.median:
            return (
                "WARN",
                f"ess_min {_fmt(new.median)} (was {_fmt(base.median)})",
            )
        return ("PASS", f"ess_min {_fmt(new.median)}")
    if q == "divergences":
        if new.median - base.median >= 5 or new.median > 10:
            return (
                "FAIL",
                f"divergences {int(new.median)} (was {int(base.median)})",
            )
        if new.median > base.median:
            return (
                "WARN",
                f"divergences {int(new.median)} (was {int(base.median)})",
            )
        return ("PASS", f"divergences {int(new.median)}")
    return ("PASS", "")


_BADGE = {"PASS": "✅", "WARN": "⚠️", "FAIL": "❌"}


def _diff_csv(name: str, base: Path, new: Path) -> str:
    out = [f"### `{name}`\n"]
    if base is None:
        out.append("_No baseline CSV — first run; recording only._\n")
        # Still emit the new rows as a single-column table for context.
        rows = _read(new)
        out.append("| section | quantity | new median [95% CrI] |")
        out.append("|---|---|---|")
        for (section, quantity), r in rows.items():
            out.append(f"| {section} | {quantity} | {_interval(r)} |")
        return "\n".join(out) + "\n"
    if new is None:
        out.append("_No new CSV — walkthrough was not built._\n")
        return "\n".join(out) + "\n"

    base_rows = _read(base)
    new_rows = _read(new)

    out.append(
        "| status | section | quantity | baseline median [95% CrI] |"
        " new median [95% CrI] | note |"
    )
    out.append("|---|---|---|---|---|---|")

    keys = sorted(set(base_rows) | set(new_rows))
    summary = {"PASS": 0, "WARN": 0, "FAIL": 0, "ONLY_BASELINE": 0, "ONLY_NEW": 0}
    for key in keys:
        section, quantity = key
        if key not in new_rows:
            summary["ONLY_BASELINE"] += 1
            out.append(
                f"| ⚠️ | {section} | {quantity} |"
                f" {_interval(base_rows[key])} | — | dropped from PR |"
            )
            continue
        if key not in base_rows:
            summary["ONLY_NEW"] += 1
            out.append(
                f"| ⚠️ | {section} | {quantity} | — |"
                f" {_interval(new_rows[key])} | new in PR |"
            )
            continue
        b = base_rows[key]
        n = new_rows[key]
        if section == "diagnostics":
            status, note = _status_diagnostic(b, n)
        else:
            status, note = _status_headline(b, n)
        summary[status] += 1
        out.append(
            f"| {_BADGE[status]} | {section} | {quantity} |"
            f" {_interval(b)} | {_interval(n)} | {note} |"
        )

    out.insert(
        1,
        "**Totals:** "
        + ", ".join(
            f"{_BADGE.get(k, '·')} {k}={v}"
            for k, v in summary.items()
            if v
        )
        + "\n",
    )
    return "\n".join(out) + "\n"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("baseline_dir", type=Path)
    p.add_argument("new_dir", type=Path)
    p.add_argument("--out", type=Path)
    args = p.parse_args()

    base_csvs = {p.name: p for p in args.baseline_dir.glob("*.csv")} \
        if args.baseline_dir.exists() else {}
    new_csvs = {p.name: p for p in args.new_dir.glob("*.csv")} \
        if args.new_dir.exists() else {}
    names = sorted(set(base_csvs) | set(new_csvs))

    if not names:
        report = "_No regression-summary CSVs found in either directory._\n"
    else:
        chunks = [
            "## 📊 Regression summary",
            "",
            "Per-walkthrough diff against the checked-in regression baseline."
            " Tolerances: 2% of the baseline median or 10% of its 95% CrI"
            " width (whichever is larger). Diagnostics use absolute"
            " thresholds (`rhat_max > 1.05`, `ess_min < 200` or"
            " > 30% drop, `divergences > 10` or > +5).",
            "",
        ]
        for name in names:
            chunks.append(_diff_csv(name, base_csvs.get(name), new_csvs.get(name)))
        report = "\n".join(chunks)

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(report)
    sys.stdout.write(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
