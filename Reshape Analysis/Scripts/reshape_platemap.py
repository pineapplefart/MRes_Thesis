# -*- coding: utf-8 -*-
"""
RESHAPE plate map loader.

Reads RESHAPE_PlateMap.xlsx and turns it into a tidy DataFrame telling the
analysis pipeline what sits in each of the 15 tray positions (A1-E3) for
every experiment.

Workbook structure expected:
    README        - human instructions, ignored here
    Experiments   - Data_Folder -> Sheet_Name lookup (header on row 3)
    Strains       - strain reference list (header on row 3)
    <one sheet per experiment> - header on row 3, EXAMPLE on row 4,
                                 15 plate rows on rows 5-19

Imperial College London
@Author: Wilson Porteus (built on Juan F Sandoval's RESHAPE pipeline)
@Date: July 2026
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pandas as pd

# Sheet layout constants — change these only if you restructure the workbook
INDEX_SHEET = "Experiments"
STRAIN_SHEET = "Strains"
HEADER_ROW = 2          # 0-based index of the header row -> Excel row 3
EXAMPLE_ROWS = 1        # number of example rows directly under the header
N_PLATES = 15

ROLE_SAMPLE = "sample"
ROLE_WATER = "water"
ROLE_CONTAM = "contamination_control"
ROLE_MINUS = "minus_variable_control"
# A position deliberately left out — a missing replicate, or a control that was not
# set up. Marking it says "this is known and intended" rather than "not filled in
# yet", so the loader stays quiet about it.
ROLE_UNUSED = "unused"
VALID_ROLES = {ROLE_SAMPLE, ROLE_WATER, ROLE_CONTAM, ROLE_MINUS, ROLE_UNUSED}

# Counts for a complete tray. A shortfall is reported for information only: not
# every tray is complete, and that is a fact about the experiment rather than an
# error in the map.
EXPECTED_COUNTS = {ROLE_SAMPLE: 12, ROLE_WATER: 1, ROLE_CONTAM: 1, ROLE_MINUS: 1}

PLATE_COLUMNS = [
    "Plate", "Grid_Row", "Grid_Col", "Role", "Strain", "Replicate",
    "Spores", "Carbon", "Nitrogen", "Include", "Baseline_Frame", "Notes",
]


def _clean(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def _yes(x) -> bool:
    return _clean(x).upper() in {"Y", "YES", "TRUE", "1"}


def _baseline(x):
    """Parse a Baseline_Frame cell. Returns an int, or None meaning 'auto'."""
    text = _clean(x).lower()
    if text in ("", "auto", "automatic", "detect"):
        return None
    try:
        value = int(float(text))
    except ValueError:
        return None
    return max(0, value)


def _offset(x) -> float:
    """Parse Start_Offset_h. Accepts 1.5, '1.5', '1:30' or '1h30'."""
    text = _clean(x)
    if text == "":
        return 0.0
    try:
        return float(text)
    except ValueError:
        pass
    match = re.match(r"^(\d+)\s*[h:]\s*(\d+)", text, re.IGNORECASE)
    if match:
        return int(match.group(1)) + int(match.group(2)) / 60.0
    return 0.0


def load_experiment_index(workbook_path) -> pd.DataFrame:
    """Return the Data_Folder -> Sheet_Name lookup table."""
    df = pd.read_excel(workbook_path, sheet_name=INDEX_SHEET, header=HEADER_ROW)
    df = df.dropna(subset=["Data_Folder", "Sheet_Name"])
    for col in ("Data_Folder", "Sheet_Name", "Carbon_Source", "Nitrogen_Source"):
        if col in df.columns:
            df[col] = df[col].map(_clean)

    df["Start_Offset_h"] = (df["Start_Offset_h"].map(_offset)
                            if "Start_Offset_h" in df.columns else 0.0)
    df["Baseline_Frame"] = (df["Baseline_Frame"].map(_baseline)
                            if "Baseline_Frame" in df.columns else None)

    # Runs sharing a Condition are biological replicates and get pooled. A blank
    # Condition means the run stands alone, so it becomes its own condition.
    if "Condition" in df.columns:
        df["Condition"] = df["Condition"].map(_clean)
    else:
        df["Condition"] = ""
    df.loc[df["Condition"] == "", "Condition"] = df["Data_Folder"]

    if "Bio_Replicate" in df.columns:
        df["Bio_Replicate"] = pd.to_numeric(df["Bio_Replicate"], errors="coerce")
    else:
        df["Bio_Replicate"] = pd.NA
    # Number unlabelled runs within each condition, in sheet order.
    missing = df["Bio_Replicate"].isna()
    if missing.any():
        df.loc[missing, "Bio_Replicate"] = (
            df[missing].groupby("Condition").cumcount() + 1
        )
    df["Bio_Replicate"] = df["Bio_Replicate"].astype(int)

    # Which component the minus-variable control leaves out. For the dextrose runs
    # that is the nitrogen source (ammonium sulfate), not the dextrose.
    if "Minus_Control_Omits" in df.columns:
        df["Minus_Control_Omits"] = df["Minus_Control_Omits"].map(_clean)
    else:
        df["Minus_Control_Omits"] = ""

    # Last usable frame, inclusive. Blank means use every frame. Set it where the
    # tail of a run is unusable — a lamp left off or the tray handled during the
    # final capture makes those frames black or blown out.
    if "Last_Frame" in df.columns:
        df["Last_Frame"] = pd.to_numeric(df["Last_Frame"], errors="coerce")
    else:
        df["Last_Frame"] = pd.NA

    if "Gelling_Agent" in df.columns:
        df["Gelling_Agent"] = df["Gelling_Agent"].map(_clean)
    else:
        df["Gelling_Agent"] = ""

    return df.reset_index(drop=True)


def load_strains(workbook_path) -> list[str]:
    df = pd.read_excel(workbook_path, sheet_name=STRAIN_SHEET, header=HEADER_ROW)
    if "Strain_ID" not in df.columns:
        return []
    return [s for s in df["Strain_ID"].map(_clean).tolist() if s]


def load_plate_map(workbook_path, strict: bool = False) -> pd.DataFrame:
    """
    Load every experiment sheet into one tidy long DataFrame.

    Returns columns:
        Experiment, Data_Folder, Sheet, Plate, Grid_Row, Grid_Col, Role,
        Strain, Replicate, Spores, Carbon, Nitrogen, Include, Notes,
        Group  (label used for averaging across replicates)

    strict=True raises on any validation problem; otherwise problems are
    printed as warnings and the offending rows are dropped.
    """
    workbook_path = Path(workbook_path)
    if not workbook_path.exists():
        raise FileNotFoundError(f"Plate map workbook not found: {workbook_path}")

    index = load_experiment_index(workbook_path)
    known_strains = set(load_strains(workbook_path))
    problems: list[str] = []   # things that look like mistakes
    notes: list[str] = []      # things that are simply true of the experiment
    frames = []

    for _, meta in index.iterrows():
        sheet = meta["Sheet_Name"]
        folder = meta["Data_Folder"]

        try:
            raw = pd.read_excel(workbook_path, sheet_name=sheet, header=HEADER_ROW)
        except ValueError:
            problems.append(f"[{folder}] sheet '{sheet}' not found in workbook")
            continue

        missing_cols = [c for c in PLATE_COLUMNS if c not in raw.columns]
        if missing_cols:
            problems.append(f"[{sheet}] missing columns: {missing_cols}")
            continue

        df = raw[PLATE_COLUMNS].iloc[EXAMPLE_ROWS:EXAMPLE_ROWS + N_PLATES].copy()
        df = df[df["Plate"].map(_clean) != ""]
        df["Plate"] = df["Plate"].map(lambda v: _clean(v).upper())
        df = df[df["Plate"] != "EXAMPLE"]

        df["Role"] = df["Role"].map(lambda v: _clean(v).lower())
        df["Strain"] = df["Strain"].map(_clean)
        df["Notes"] = df["Notes"].map(_clean)
        df["Spores"] = df["Spores"].map(_yes)
        df["Include"] = df["Include"].map(_yes)
        df["Replicate"] = pd.to_numeric(df["Replicate"], errors="coerce")
        df["Carbon"] = df["Carbon"].map(_clean)
        df["Nitrogen"] = df["Nitrogen"].map(_clean)

        # Per-plate baseline overrides the experiment default; None means auto-detect.
        sheet_baseline = meta.get("Baseline_Frame")
        df["Baseline_Frame"] = pd.to_numeric(
            df["Baseline_Frame"].map(_baseline), errors="coerce"
        ).astype("float64")
        if sheet_baseline is not None:
            df["Baseline_Frame"] = df["Baseline_Frame"].fillna(float(sheet_baseline))
        df["Start_Offset_h"] = float(meta.get("Start_Offset_h") or 0.0)

        # --- validation -------------------------------------------------
        blank = df[df["Role"] == ""]["Plate"].tolist()
        if len(blank) == len(df):
            # Sheet not started yet — say so once and move on quietly.
            problems.append(f"[{sheet}] not filled in yet — skipped")
            continue
        if blank:
            problems.append(f"[{sheet}] Role not filled in for plates: {blank}")
        df = df[df["Role"] != ""]

        bad_role = df[~df["Role"].isin(VALID_ROLES)]
        if len(bad_role):
            problems.append(
                f"[{sheet}] unrecognised Role values: "
                f"{sorted(set(bad_role['Role']))} on plates {bad_role['Plate'].tolist()}"
            )
            df = df[df["Role"].isin(VALID_ROLES)]

        counts = df["Role"].value_counts().to_dict()
        shortfalls = []
        for role, want in EXPECTED_COUNTS.items():
            got = counts.get(role, 0)
            if got != want:
                shortfalls.append(f"{got}/{want} {role}")
        if shortfalls:
            notes.append(f"[{sheet}] incomplete tray: " + ", ".join(shortfalls))

        samples = df[df["Role"] == ROLE_SAMPLE]
        no_strain = samples[samples["Strain"] == ""]["Plate"].tolist()
        if no_strain:
            problems.append(f"[{sheet}] sample plates with no Strain: {no_strain}")

        no_rep = samples[samples["Replicate"].isna()]["Plate"].tolist()
        if no_rep:
            problems.append(f"[{sheet}] sample plates with no Replicate: {no_rep}")

        if known_strains:
            unknown = sorted(set(samples["Strain"]) - known_strains - {""})
            if unknown:
                problems.append(f"[{sheet}] strains not listed on the Strains sheet: {unknown}")

        dupes = (
            samples[samples["Strain"] != ""]
            .groupby(["Strain", "Replicate"], dropna=False)
            .size()
        )
        dupes = dupes[dupes > 1]
        if len(dupes):
            problems.append(f"[{sheet}] duplicate Strain/Replicate pairs: {dupes.to_dict()}")

        dup_plates = df[df["Plate"].duplicated()]["Plate"].tolist()
        if dup_plates:
            problems.append(f"[{sheet}] duplicate plate positions: {dup_plates}")

        # --- labels -----------------------------------------------------
        df.insert(0, "Sheet", sheet)
        df.insert(0, "Data_Folder", folder)
        df.insert(0, "Experiment", folder)
        df["Condition"] = meta.get("Condition") or folder
        df["Bio_Replicate"] = int(meta.get("Bio_Replicate") or 1)
        df["Minus_Control_Omits"] = meta.get("Minus_Control_Omits") or ""
        last = meta.get("Last_Frame")
        df["Last_Frame"] = float(last) if pd.notna(last) else float("nan")
        df["Gelling_Agent"] = meta.get("Gelling_Agent") or ""
        df["Group"] = df.apply(_group_label, axis=1)
        frames.append(df)

    if notes:
        print("Plate map notes (not errors):\n  - " + "\n  - ".join(notes))

    if problems:
        message = "Plate map validation:\n  - " + "\n  - ".join(problems)
        if strict:
            raise ValueError(message)
        print(message)

    if not frames:
        raise ValueError(
            "No usable plate rows found. Fill in the Role column on at least one "
            f"experiment sheet of {workbook_path.name}."
        )

    out = pd.concat(frames, ignore_index=True)
    return out


def _group_label(row) -> str:
    """
    Label used to average replicates together. Keyed on Condition, not Experiment,
    so separate runs of the same condition pool into one group.
    """
    condition = row.get("Condition") or row["Experiment"]
    role = row["Role"]
    if role == ROLE_SAMPLE:
        return f"{condition} | {row['Strain']}"
    if role == ROLE_MINUS:
        # Name the component that was left out — "minus_variable_control" alone does
        # not say what the plate is missing, and it differs between runs.
        omits = row.get("Minus_Control_Omits") or ""
        if omits:
            return f"{condition} | no {omits.lower()}"
    return f"{condition} | {role}"


def plate_lookup(plate_map: pd.DataFrame, experiment: str) -> dict:
    """dict of plate code -> metadata row, for one experiment."""
    sub = plate_map[plate_map["Experiment"] == experiment]
    return {r["Plate"]: r.to_dict() for _, r in sub.iterrows()}


def resolve_experiment(plate_map: pd.DataFrame, folder_name: str) -> str | None:
    """Match a data folder (or zip stem) on disk to an Experiment in the map."""
    name = folder_name.strip()
    exact = plate_map[plate_map["Data_Folder"].str.lower() == name.lower()]
    if len(exact):
        return exact.iloc[0]["Experiment"]
    return None


if __name__ == "__main__":
    try:
        BASE_DIR = Path(__file__).resolve().parent.parent
    except NameError:
        BASE_DIR = Path(os.getcwd())

    wb = BASE_DIR / "RESHAPE_PlateMap.xlsx"
    pm = load_plate_map(wb)
    print(f"\nLoaded {len(pm)} plate rows across {pm['Experiment'].nunique()} experiment(s).")
    print(pm.groupby(["Experiment", "Role"]).size())
