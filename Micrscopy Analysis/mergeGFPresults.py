"""
Merge GFP Quantification Results
=================================
Scans a root folder for the per-strain/condition summary CSVs produced by
the Fiji macro (GFP_quantification_macro.ijm), e.g.:

    2.1_agar_summary.csv
    3.1_lc_t9.5_summary.csv
    2.1_tapelift_summary.csv

...and merges them into a single master table with Strain, Condition,
Timepoint (parsed out of "lc_tX" conditions, blank otherwise), and all the
original Analyze Particles summary columns (Count, Total Area, %Area, Mean).

Usage:
    python merge_gfp_results.py /path/to/root_folder

Output:
    <root_folder>/merged_gfp_results.csv
    <root_folder>/merged_gfp_results.xlsx
"""

import sys
import re
from pathlib import Path

import pandas as pd

# Known strains and non-timepoint condition names -- edit here if your
# naming changes.
STRAINS = ["1.1", "2.1", "3.1"]
NON_TIMEPOINT_CONDITIONS = {"agar", "glass", "glass_pbs", "plug", "tapelift"}


def parse_filename(filename: str):
    """
    Given a filename like '2.1_lc_t9.5_summary.csv' or '3.1_agar_summary.csv',
    return (strain, condition, timepoint_or_None).
    """
    name = filename.replace("_summary.csv", "")

    strain = None
    for s in STRAINS:
        if name.startswith(s + "_"):
            strain = s
            condition = name[len(s) + 1:]
            break
    else:
        return None, name, None

    timepoint = None
    if condition in NON_TIMEPOINT_CONDITIONS:
        pass  # no timepoint for these
    else:
        # expects condition like "lc_t0", "lc_t3", "lc_t9.5", "lc_t24"
        match = re.match(r"lc_t([\d.]+)", condition)
        if match:
            timepoint = match.group(1)

    return strain, condition, timepoint


def main(root_folder: str):
    root = Path(root_folder)
    if not root.is_dir():
        print(f"Error: {root_folder} is not a valid directory.")
        sys.exit(1)

    csv_files = sorted(root.glob("*_summary.csv"))
    if not csv_files:
        print(f"No '*_summary.csv' files found in {root_folder}.")
        print("Make sure you've run the Fiji macro first.")
        sys.exit(1)

    all_rows = []
    for csv_path in csv_files:
        strain, condition, timepoint = parse_filename(csv_path.name)
        if strain is None:
            print(f"  Skipping (couldn't parse strain): {csv_path.name}")
            continue

        df = pd.read_csv(csv_path)
        df.insert(0, "Strain", strain)
        df.insert(1, "Condition", condition)
        df.insert(2, "Timepoint", timepoint)
        df.insert(3, "SourceFile", csv_path.name)
        all_rows.append(df)
        print(f"  Loaded {csv_path.name}: {len(df)} row(s) -> strain={strain}, condition={condition}, timepoint={timepoint}")

    if not all_rows:
        print("No valid CSVs were parsed. Nothing to merge.")
        sys.exit(1)

    merged = pd.concat(all_rows, ignore_index=True)

    out_csv = root / "merged_gfp_results.csv"
    out_xlsx = root / "merged_gfp_results.xlsx"
    merged.to_csv(out_csv, index=False)
    merged.to_excel(out_xlsx, index=False)

    print(f"\nMerged {len(csv_files)} files into {len(merged)} total rows.")
    print(f"Saved: {out_csv}")
    print(f"Saved: {out_xlsx}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python merge_gfp_results.py /path/to/root_folder")
        sys.exit(1)
    main(sys.argv[1])