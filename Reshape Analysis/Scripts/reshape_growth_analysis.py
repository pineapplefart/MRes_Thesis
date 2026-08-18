# -*- coding: utf-8 -*-
"""
RESHAPE radial growth analysis — image-sequence version.

Adapted from Juan F Sandoval's RESHAPE video pipeline for RI92 timelapse data
that is exported as one JPG sequence per plate rather than one video per plate.

Expected data layout (either a folder or a .zip, both work):

    Data/
      Dextrose Urea.zip                 <- or an extracted folder of the same name
        Dextrose Urea/
          A1/ RIS1_0_TL_10-PlateA1-0.jpg ... -95.jpg
          A2/ ...
          ...
          E3/ ...

The 15 positions A1-E3 are mapped to experimental conditions by
RESHAPE_PlateMap.xlsx (see reshape_platemap.py).

What it does per plate:
  1. Detect the agar circle once, on a reference frame (Hough, as in Juan's code).
  2. Mask + crop every frame to that circle.
  3. Background-subtract against frame 0 to segment mycelium.
  4. Record colony area, equivalent radius and mean RGB/luminance per frame.
  5. Average replicates by strain and plot growth curves, with controls handled
     separately.

Frames are processed in memory — no intermediate PNGs are written unless
SAVE_QC_FRAMES is on — so a 23 GB dataset does not become a 200 GB one.

Imperial College London
@Author: Wilson Porteus (built on Juan F Sandoval's RESHAPE pipeline)
@Date: July 2026
"""

from __future__ import annotations

import hashlib
import io
import math
import os
import pickle
import re
import zipfile
from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from reshape_summary import (
    kr_matrix,
    plot_kr_vs_mu,
    plot_minus_controls,
    minus_control_check,
    per_plate_rates,
    plot_kr_heatmap,
    strain_substrate_summary,
    write_summary_workbook,
)
from reshape_platemap import (
    ROLE_CONTAM,
    ROLE_UNUSED,
    ROLE_MINUS,
    ROLE_SAMPLE,
    ROLE_WATER,
    load_plate_map,
    plate_lookup,
    resolve_experiment,
)

# ---------------------------------------------------------------- paths
try:
    BASE_DIR = Path(__file__).resolve().parent.parent
except NameError:
    BASE_DIR = Path(os.getcwd())

INPUT_DIR = BASE_DIR / "Data"
PLATEMAP_PATH = BASE_DIR / "RESHAPE_PlateMap.xlsx"
OUTPUT_DIR = BASE_DIR / "Outputs"
OUTPUT_DIR.mkdir(exist_ok=True)
CACHE_DIR = OUTPUT_DIR / "cache"
CACHE_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------- settings
# Restrict the run to specific experiments (folder names). Empty = all of them.
ONLY_EXPERIMENTS: list[str] = []

# Roles skipped entirely: not processed, not plotted, not in any CSV. Water plates
# carry no medium and so have no growth to measure. Set to () to bring them back;
# to drop an individual plate instead, set Include = N on its row in the workbook.
# "unused" marks a tray position that was deliberately not set up.
EXCLUDE_ROLES: tuple[str, ...] = (ROLE_WATER, ROLE_UNUSED)

FRAME_INTERVAL_HOURS = 0.5     # RI92 captured every 30 min
PLATE_DIAMETER_MM = 90.0
REFERENCE_FRAME_IDX = 2        # frame used to fit the agar circle

# --- condensation / baseline frame -----------------------------------------
# Growth is measured against a "clean" background frame. On plates that start
# with condensation, frame 0 is not clean: the condensation clearing shows up as
# a large apparent colony. When a plate's Baseline_Frame cell is blank the
# pipeline finds the first frame at which the image has stopped changing rapidly
# and uses that. Frames before the baseline are dropped from the series.
AUTO_BASELINE = True
AUTO_BASELINE_MAX_FRAME = 8      # never look past 4 h for a settled frame
AUTO_BASELINE_REF_FRAME = 20     # settled reference: 10 h, after condensation, before growth
AUTO_BASELINE_TOLERANCE = 1.5    # settled = within this multiple of the plateau level
AUTO_BASELINE_MIN_RATIO = 2.0    # frame 0 must be this much worse than the plateau
                                 # before we accept that condensation was present

# Source frames are 3040x3040. Analysing at full resolution is ~16x slower and
# changes area fractions by well under a percent. 0.25 -> 760x760.
DOWNSCALE = 0.25

# Agar detection (Hough) — radii are fractions of the frame's short side
HOUGH_MIN_RADIUS_FRAC = 0.35
HOUGH_MAX_RADIUS_FRAC = 0.48
HOUGH_DP = 1.2
HOUGH_MIN_DIST = 200
HOUGH_PARAM1 = 100
HOUGH_PARAM2 = 100
FALLBACK_RADIUS_FRAC = 0.45    # used if Hough finds nothing

# --- mycelium segmentation --------------------------------------------------
# Mycelium is brighter than the medium on every substrate here, so we threshold a
# signed rise in CIE lightness against the baseline frame rather than any colour
# change. Three things make that hold up across substrates:
#
#   RIM_KEEP          the Hough circle lands on the OUTER dish rim, and the wall
#                     plus meniscus occupy the outer ~20% of that radius. Measured
#                     on water controls, the artifact sits in the 0.8-0.9 band, so
#                     everything beyond RIM_KEEP is discarded. Without this the
#                     rim reflection is picked up as a colony.
#   drift correction  opaque media (cellulose, molasses) darken markedly over 48 h.
#                     Each frame the medium's own shift is measured from the pixels
#                     not yet claimed by the colony and subtracted, so a global
#                     change cannot flood the plate.
#   drift freezing    once the colony leaves less than MIN_MEDIUM_FRAC of visible
#                     medium there is nothing left to measure, and re-estimating
#                     would lock onto the colony itself and erase it. The last
#                     reliable value is held instead.
SEG_VERSION = 7                  # bump to invalidate cached results
RIM_KEEP = 0.78                  # fraction of detected radius analysed
MIN_DELTA_L = 3.0                # floor on brightness rise, 8-bit L units
MAD_K = 2.5                      # threshold = drift + MAD_K * robust sigma
MONOTONIC = True                 # a colony on agar does not shrink

# Segmentation runs in TWO PASSES over each plate.
#
# Pass 1 sweeps forward, estimating the medium drift and thresholding each frame,
# and unions the results into the colony's full footprint. Pass 2 replays the frames
# restricted to that footprint and accumulates, so the mask can only grow.
#
# A single streaming pass cannot do this. Deciding frame by frame which blob is the
# colony lets the choice jump between objects, which makes the mask shrink — the
# negative growth rates. Locking onto the first blob instead fixes the shrinking but
# freezes the plate on whatever it latched onto first: Molasses A1 locked onto an
# 8950 px artifact at frame 10 and stayed there for the remaining 43 hours while the
# real colony, a separate component, was never picked up. Size cannot separate the
# two — that artifact was larger than many genuine early colonies.
#
# Anchoring on the footprint works because the colony is unambiguous in the frames
# where it is largest, and growth on agar is spatially nested: mask(t) is contained
# in mask(T).
FOOTPRINT_DILATE_FRAC = 0.01     # slack around the footprint, fraction of plate radius
SEED_CORE_FRAC = 0.35            # inoculum region = this deep into the footprint

# Freezing the drift estimate. On a plate the mycelium eventually covers, the
# "medium pool" is really colony, and re-estimating from it walks the drift upward
# until the threshold exceeds the colony's own brightness and growth stops being
# detected — Molasses A1 drifted to +14 L and stalled at 58% on a fully covered
# plate. The estimate is therefore frozen while the pool is still mostly medium.
#
# Tempting alternative, do not use: capping the drift so it cannot rise. Some media
# (the blue agars) legitimately brighten, and capping makes the whole plate read as
# colony — it put a 43% mask on a visibly clean Dextrose Urea plate at 24 h.
MIN_MEDIUM_FRAC = 0.60
BLUR_KERNEL = 5
MORPH_KERNEL = 5
MIN_COMPONENT_AREA_FRAC = 2e-5   # fraction of agar area, scale-independent

# Loose aerial hyphae threshold as disconnected filaments with gaps between them.
# Each plate carries exactly one colony, so the mask is closed up and returned as a
# single filled region.
#   GAP_CLOSE_FRAC  morphological closing radius, as a fraction of the plate radius.
#                   Bridges gaps narrower than this between filaments of the same
#                   colony. Closing is dilate-then-erode, so the outer boundary is
#                   largely preserved; it fills gaps rather than growing the colony.
#   FILL_HOLES      fill anything enclosed by the colony outline. The interior of a
#                   colony is colonised even where individual hyphae are sparse.
GAP_CLOSE_FRAC = 0.06
FILL_HOLES = True

# Spores are pipetted into the middle, so the colony starts at the centre and grows
# outward. Components are therefore scored on how close their PIXELS get to the
# centre, not on their centroid: a long arc hugging the rim has a near-central
# centroid, which is how rim reflections used to win. Radii below are fractions of
# the plate radius.
NEAR_WEIGHT = 0.8                # penalty on the component's closest approach
MEAN_WEIGHT = 0.4                # penalty on its average distance from centre
MAX_START_RADIUS = 0.55          # reject anything never reaching within this of centre

# Output control
USE_CACHE = True               # reuse per-experiment pickles from a previous run
SAVE_QC_FRAMES = True          # write a few overlay PNGs per plate for eyeballing
QC_FRAME_INDICES = (0, 24, 48, 71, 95)
PLOTS = True
EXPORT_CSV = True

plt.rc("font", family="sans-serif", weight="normal", size=13)
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

FRAME_RE = re.compile(r"-Plate([A-Z]\d+)-(\d+)\.jpe?g$", re.IGNORECASE)


# ================================================================ sources
class PlateSource:
    """Uniform access to plate frames, whether they live in a folder or a zip."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.is_zip = self.path.suffix.lower() == ".zip"
        self._zf = None
        self._index: dict[str, list[tuple[int, str]]] = {}
        self._build_index()

    @property
    def name(self) -> str:
        return self.path.stem if self.is_zip else self.path.name

    def _build_index(self):
        if self.is_zip:
            self._zf = zipfile.ZipFile(self.path)
            names = self._zf.namelist()
        else:
            names = [str(p.relative_to(self.path)) for p in self.path.rglob("*.jp*g")]

        for n in names:
            if n.endswith("/"):
                continue
            m = FRAME_RE.search(n)
            if not m:
                continue
            plate = m.group(1).upper()
            idx = int(m.group(2))
            self._index.setdefault(plate, []).append((idx, n))

        for plate in self._index:
            self._index[plate].sort(key=lambda t: t[0])

    def plates(self) -> list[str]:
        return sorted(self._index.keys())

    def frame_indices(self, plate: str) -> list[int]:
        return [i for i, _ in self._index[plate]]

    def read_frame(self, plate: str, position: int) -> np.ndarray:
        _idx, member = self._index[plate][position]
        if self.is_zip:
            buf = np.frombuffer(self._zf.read(member), dtype=np.uint8)
            img = cv2.imdecode(buf, cv2.IMREAD_COLOR)
        else:
            img = cv2.imread(str(self.path / member), cv2.IMREAD_COLOR)
        if img is None:
            raise ValueError(f"Could not decode frame: {member}")
        if DOWNSCALE and DOWNSCALE != 1.0:
            img = cv2.resize(img, None, fx=DOWNSCALE, fy=DOWNSCALE,
                             interpolation=cv2.INTER_AREA)
        return img

    def close(self):
        if self._zf is not None:
            self._zf.close()


def discover_sources(input_dir: Path) -> list[PlateSource]:
    """Prefer an extracted folder over a zip of the same name."""
    folders = {p.name: p for p in input_dir.iterdir()
               if p.is_dir() and not p.name.startswith(".")}
    zips = {p.stem: p for p in input_dir.glob("*.zip")}

    chosen: dict[str, Path] = dict(zips)
    chosen.update(folders)          # folder wins if both exist

    sources = []
    for name in sorted(chosen):
        try:
            src = PlateSource(chosen[name])
        except zipfile.BadZipFile:
            print(f"[SKIP] {name}: not a readable zip")
            continue
        if src.plates():
            sources.append(src)
        else:
            src.close()
    return sources


# ================================================================ detection
def detect_agar_circle(image_bgr: np.ndarray) -> tuple[int, int, int]:
    """Fit the agar circle. Falls back to a centred circle if Hough finds none."""
    h, w = image_bgr.shape[:2]
    short = min(h, w)
    min_r = int(short * HOUGH_MIN_RADIUS_FRAC)
    max_r = int(short * HOUGH_MAX_RADIUS_FRAC)

    gray = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2GRAY)
    gray = cv2.GaussianBlur(gray, (9, 9), 0)

    circles = cv2.HoughCircles(
        gray, cv2.HOUGH_GRADIENT,
        dp=HOUGH_DP,
        minDist=max(1, int(HOUGH_MIN_DIST * (DOWNSCALE or 1.0))),
        param1=HOUGH_PARAM1,
        param2=HOUGH_PARAM2,
        minRadius=min_r,
        maxRadius=max_r,
    )

    img_center = np.array([w / 2.0, h / 2.0])

    if circles is None:
        r = int(short * FALLBACK_RADIUS_FRAC)
        return int(img_center[0]), int(img_center[1]), r

    circles = np.round(circles[0]).astype(int)
    best, best_score = None, None
    for x, y, r in circles:
        dist = np.linalg.norm(np.array([x, y]) - img_center)
        score = dist - 0.2 * r            # near-centre and large is good
        if best_score is None or score < best_score:
            best_score, best = score, (int(x), int(y), int(r))
    return best


def crop_to_circle(frame_bgr: np.ndarray, cx: int, cy: int, r: int):
    """Return (bgr_crop, valid_mask) for the disc, both cropped to its bbox."""
    h, w = frame_bgr.shape[:2]
    x1, y1 = max(cx - r, 0), max(cy - r, 0)
    x2, y2 = min(cx + r, w), min(cy + r, h)

    crop = frame_bgr[y1:y2, x1:x2]
    mask = np.zeros((y2 - y1, x2 - x1), dtype=np.uint8)
    cv2.circle(mask, (cx - x1, cy - y1), r, 255, thickness=-1)
    return crop, mask > 0


def radius_map(shape) -> np.ndarray:
    """Distance from the crop centre, as a fraction of the detected plate radius."""
    h, w = shape
    ys, xs = np.mgrid[0:h, 0:w]
    rr = np.sqrt((xs - (w - 1) / 2.0) ** 2 + (ys - (h - 1) / 2.0) ** 2)
    return rr / (min(h, w) / 2.0)


def agar_mask(shape, valid: np.ndarray, keep: float = RIM_KEEP) -> np.ndarray:
    """The analysed region: the disc with the dish wall and meniscus removed."""
    return valid & (radius_map(shape) <= keep)


def mean_colour(crop_bgr: np.ndarray, valid: np.ndarray):
    b, g, r = cv2.split(crop_bgr)
    mean_r = float(np.mean(r[valid]))
    mean_g = float(np.mean(g[valid]))
    mean_b = float(np.mean(b[valid]))
    lum = 0.2126 * mean_r + 0.7152 * mean_g + 0.0722 * mean_b
    return mean_r, mean_g, mean_b, float(lum)


def area_metrics(colony_mask: np.ndarray, agar: np.ndarray, mm_per_px: float):
    """
    Areas in mm2 are calibrated from the detected dish rim (PLATE_DIAMETER_MM across),
    not from the pixel count of the analysed region, so narrowing the analysed disc
    to exclude the wall does not change the physical scale.
    """
    colony_px = int(np.sum(colony_mask))
    agar_px = int(np.sum(agar))
    if agar_px == 0:
        raise ValueError("Agar mask contains zero valid pixels.")

    colony_area_mm2 = colony_px * (mm_per_px ** 2)

    return {
        "colony_pixels": colony_px,
        "agar_pixels": agar_px,
        "area_percent": 100.0 * colony_px / agar_px,
        "area_mm2": colony_area_mm2,
        "equivalent_radius_px": math.sqrt(colony_px / math.pi) if colony_px else 0.0,
        "equivalent_radius_mm": math.sqrt(colony_area_mm2 / math.pi) if colony_area_mm2 else 0.0,
    }


def _lab_lightness(bgr: np.ndarray) -> np.ndarray:
    k = BLUR_KERNEL + (BLUR_KERNEL % 2 == 0)
    blurred = cv2.GaussianBlur(bgr, (k, k), 0)
    return cv2.cvtColor(blurred, cv2.COLOR_BGR2LAB)[:, :, 0].astype(np.float32)


def _robust_centre_scale(values: np.ndarray, rounds: int = 3):
    """Median and MAD-sigma, iteratively trimming the bright tail (young colony)."""
    v = values
    centre = scale = 0.0
    for _ in range(rounds):
        if v.size < 200:
            break
        centre = float(np.median(v))
        scale = max(float(np.median(np.abs(v - centre))) * 1.4826, 0.5)
        v = v[v <= centre + MAD_K * scale]
    return centre, max(scale, 0.5)


class PlateSegmenter:
    """
    Two-pass segmentation of one plate's frames. See the notes beside
    FOOTPRINT_DILATE_FRAC for why a single streaming pass does not work.

    Usage:
        seg = PlateSegmenter(base_crop, valid)
        for crop in crops:
            seg.add(crop)          # stores lightness only, not the full frame
        masks, diagnostics = seg.run()
    """

    def __init__(self, base_crop: np.ndarray, valid: np.ndarray):
        self.agar = agar_mask(base_crop.shape[:2], valid)
        self.L_ref = _lab_lightness(base_crop)
        self._L: list[np.ndarray] = []

        self._morph = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (MORPH_KERNEL, MORPH_KERNEL))
        self._halo = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (15, 15))
        self._rmap = radius_map(base_crop.shape[:2])

        radius_px = min(base_crop.shape[:2]) / 2.0
        gap = max(3, int(round(GAP_CLOSE_FRAC * radius_px)) | 1)
        self._gap = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (gap, gap))
        grow = max(3, int(round(FOOTPRINT_DILATE_FRAC * radius_px)) | 1)
        self._foot = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (grow, grow))

        self.frozen_at = None
        self.footprint = None
        self.seed = None

    # -- frame intake ----------------------------------------------------
    def add(self, crop_bgr: np.ndarray):
        """Store only the lightness plane; keeping full frames would be ~4x the RAM."""
        self._L.append(_lab_lightness(crop_bgr).astype(np.uint8))

    # -- helpers ---------------------------------------------------------
    def _clean(self, raw_bool):
        u8 = (raw_bool.astype(np.uint8) * 255)
        out = cv2.morphologyEx(u8, cv2.MORPH_OPEN, self._morph)
        return cv2.morphologyEx(out, cv2.MORPH_CLOSE, self._gap)

    def _single_colony(self, binary_u8):
        """Pick the one component most consistent with growth from the inoculum."""
        n, labels, stats, _ = cv2.connectedComponentsWithStats(binary_u8, connectivity=8)
        kept = np.zeros_like(binary_u8)
        if n <= 1:
            return kept
        min_area = max(20, int(MIN_COMPONENT_AREA_FRAC * self.agar.sum()))
        best, best_score = None, None
        for label in range(1, n):
            area = stats[label, cv2.CC_STAT_AREA]
            if area < min_area:
                continue
            radii = self._rmap[labels == label]
            nearest = float(radii.min())
            average = float(radii.mean())
            if nearest > MAX_START_RADIUS:
                continue
            score = (math.sqrt(area)
                     * (1.0 - NEAR_WEIGHT * min(nearest, 1.0))
                     * (1.0 - MEAN_WEIGHT * min(average, 1.0)))
            if best_score is None or score > best_score:
                best_score, best = score, label
        if best is not None:
            kept[labels == best] = 255
        return kept

    @staticmethod
    def _component_containing(mask_bool, seed_bool):
        """The single connected component of `mask` that overlaps `seed`."""
        u8 = (mask_bool.astype(np.uint8) * 255)
        n, labels = cv2.connectedComponents(u8, connectivity=8)
        out = np.zeros_like(u8)
        if n <= 1:
            return out
        hits = labels[seed_bool]
        hits = hits[hits > 0]
        if hits.size == 0:
            return out
        out[labels == np.bincount(hits).argmax()] = 255
        return out

    @staticmethod
    def _fill_holes(binary_u8):
        contours, _ = cv2.findContours(binary_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        filled = np.zeros_like(binary_u8)
        cv2.drawContours(filled, contours, -1, 255, thickness=cv2.FILLED)
        return filled

    # -- the two passes --------------------------------------------------
    def _threshold_pass(self):
        """
        Sweep forward, tracking the medium drift, and return one raw mask per frame.

        The drift is measured from the pixels not yet claimed by anything detected so
        far, and frozen once too little medium remains to measure — otherwise the
        estimate is taken from the colony itself and walks upward until the colony
        stops being detected.
        """
        raws, picks, diags = [], [], []
        claimed = np.zeros(self.agar.shape, dtype=bool)
        drift, sigma = 0.0, 1.0
        have_estimate = False

        for idx, L in enumerate(self._L):
            delta = L.astype(np.float32) - self.L_ref
            pool = self.agar & ~(cv2.dilate(claimed.astype(np.uint8), self._halo) > 0)
            pool_frac = pool.sum() / max(1, self.agar.sum())

            if pool_frac >= MIN_MEDIUM_FRAC:
                drift, sigma = _robust_centre_scale(delta[pool])
                have_estimate = True
            elif not have_estimate:
                drift, sigma = _robust_centre_scale(delta[self.agar])
            elif self.frozen_at is None:
                self.frozen_at = idx

            threshold = max(drift + MAD_K * sigma, drift + MIN_DELTA_L)
            raw = (delta >= threshold) & self.agar
            picked = self._single_colony(self._clean(raw)) > 0

            raws.append(raw)
            picks.append(picked)
            # Claim only the selected colony, never every thresholded pixel. Claiming
            # raw detections lets scattered noise eat the medium pool on a plate with
            # no colony, which freezes the drift early on a bad value and collapses
            # the threshold — blank controls came out at 78% of the plate.
            claimed |= picked
            diags.append({"drift": drift, "sigma": sigma,
                          "threshold": threshold, "pool_frac": pool_frac})
        return raws, picks, diags

    def _find_footprint(self, picks):
        """
        The colony's full extent, taken from the SINGLE frame where its selected
        component is largest — normally near the end, where the colony is clearest.

        Not the union over all frames: on a plate with no colony that unions 96
        frames of noise into a large blob, and everything inside it then accumulates.
        A blank cellulose control read 79% of the plate that way. One frame's
        component cannot accumulate noise, so a blank plate yields a blank footprint.
        """
        best_mask, best_area = None, 0
        for picked in picks:
            area = int(picked.sum())
            if area > best_area:
                best_area, best_mask = area, picked
        if best_mask is None or best_area == 0:
            return np.zeros(self.agar.shape, dtype=bool)
        grown = cv2.dilate((best_mask.astype(np.uint8) * 255), self._foot) > 0
        footprint = grown & self.agar

        # Inoculum region: the core of the footprint, taken as the pixels at least
        # SEED_CORE_FRAC of the way to its deepest point. A whole region rather than
        # a single point, because sparse early hyphae can straddle a gap and miss a
        # point seed entirely; and taken from the footprint rather than from the
        # first detection, because the first detection is sometimes a speck the
        # colony has not reached yet, which freezes the plate on it.
        dist = cv2.distanceTransform((footprint.astype(np.uint8) * 255), cv2.DIST_L2, 5)
        if dist.max() > 0:
            self.seed = dist >= SEED_CORE_FRAC * dist.max()
        return footprint

    def run(self):
        """Return (masks, diagnostics), one per added frame."""
        raws, picks, diags = self._threshold_pass()
        self.footprint = self._find_footprint(picks)

        masks = []
        accumulated = np.zeros(self.agar.shape, dtype=bool)
        for raw, diag in zip(raws, diags):
            inside = raw & self.footprint
            # Accumulate only what survives cleaning and component selection, so
            # stray specks never enter the running union.
            picked = self._single_colony(self._clean(inside)) > 0

            if MONOTONIC:
                accumulated |= picked
                current = accumulated
            else:
                current = picked
            # Close first, so filaments of the same colony that were detected on
            # different frames read as one region rather than as islands.
            merged = cv2.morphologyEx((current.astype(np.uint8) * 255),
                                      cv2.MORPH_CLOSE, self._gap) > 0

            # Keep the component containing the inoculum. All the spores went in as
            # one droplet, so there is exactly one colony and it covers that region
            # throughout.
            #
            # This is also what makes the mask monotonic. If A is contained in B,
            # the component of A containing the seed is contained in the component
            # of B containing it — so as `accumulated` grows, so does the reported
            # colony, and it can never split in two. Unioning with the previous
            # frame instead (the earlier approach) guaranteed monotonicity but could
            # weld on a second, disconnected region — the two-mask problem.
            # Anchor on the previous mask once there is one, falling back to the
            # inoculum core on the first frames. Anchoring on the core throughout is
            # not enough: _component_containing keeps whichever component overlaps
            # the anchor most, which can hand the colony over to a newly merged
            # region and drop the previous mask. The previous mask is connected and
            # wholly inside `merged`, so the component holding it is a superset —
            # growth guaranteed, one region guaranteed.
            anchor = masks[-1] if masks and masks[-1].any() else self.seed
            grown = (self._component_containing(merged, anchor)
                     if anchor is not None
                     else np.zeros(self.agar.shape, dtype=np.uint8))
            mask = self._fill_holes(grown) > 0 if FILL_HOLES else grown > 0

            if MONOTONIC:
                # Feed the finished mask back into the accumulator rather than
                # unioning it with the previous frame's mask. Hole filling is not
                # quite monotone on its own (a hole can open onto the outside as the
                # colony grows, un-filling it), and unioning the two masks can leave
                # two disconnected regions — they may touch the seed at different
                # points. Feeding back keeps one connected region and guarantees
                # growth, because next frame's component is taken from a set that
                # already contains this one.
                accumulated |= mask

            masks.append(mask)
            diag["footprint_frac"] = float(self.footprint.sum()) / max(1, self.agar.sum())
        return masks, diags


def qc_overlay(crop_bgr, colony_mask, alpha=0.45):
    out = crop_bgr.copy()
    red = out.copy()
    red[colony_mask] = (0, 0, 255)
    out = cv2.addWeighted(red, alpha, out, 1 - alpha, 0)
    contours, _ = cv2.findContours((colony_mask.astype(np.uint8) * 255),
                                   cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    cv2.drawContours(out, contours, -1, (0, 0, 255), 2)
    return out


# ================================================================ baseline
def find_settled_frame(source: PlateSource, plate: str, cx: int, cy: int, r: int) -> tuple[int, dict]:
    """
    Find the first frame at which condensation has cleared.

    Each early frame is compared against a later reference frame that is past any
    condensation but still before growth starts (frame 20 = 10 h by default).
    A plate with condensation shows a large difference at frame 0 that decays and
    then flattens off; a clean plate is flat from the start. We take the plateau
    level to be the smallest difference in the search range and accept the first
    frame that comes within AUTO_BASELINE_TOLERANCE of it.

    If growth has already begun by the reference frame, every early frame differs
    from it by a similar amount, the plateau equals the frame-0 difference and
    this returns 0 — the safe answer.

    Returns (frame index, diagnostics).
    """
    n = len(source.frame_indices(plate))
    look = min(AUTO_BASELINE_MAX_FRAME, n - 2)
    ref_pos = min(AUTO_BASELINE_REF_FRAME, n - 1)
    if look < 1 or ref_pos <= look:
        return 0, {"reason": "sequence too short to detect condensation"}

    def disc_gray(pos: int) -> np.ndarray:
        crop, valid = crop_to_circle(source.read_frame(plate, pos), cx, cy, r)
        g = cv2.cvtColor(crop, cv2.COLOR_BGR2GRAY).astype(np.float32)
        g[~valid] = np.nan
        return g

    reference = disc_gray(ref_pos)
    diffs = {}
    for pos in range(0, look + 1):
        diffs[pos] = float(np.nanmean(np.abs(disc_gray(pos) - reference)))

    plateau = max(min(diffs.values()), 1e-6)
    threshold = AUTO_BASELINE_TOLERANCE * plateau

    # A plate that merely drifts slowly (illumination, thermal) also decays, but
    # never starts far above its plateau. Only treat frame 0 as unusable when it
    # is clearly elevated, otherwise nothing is gained by discarding frames.
    if diffs[0] < AUTO_BASELINE_MIN_RATIO * plateau:
        chosen = 0
    else:
        chosen = next((p for p in range(0, look + 1) if diffs[p] <= threshold), 0)

    diag = {
        "reference_frame": ref_pos,
        "plateau": round(plateau, 3),
        "threshold": round(threshold, 3),
        "start_ratio": round(diffs[0] / plateau, 2),
        "diffs": {k: round(v, 2) for k, v in diffs.items()},
    }
    return chosen, diag


# ================================================================ per plate
def process_plate(source: PlateSource, plate: str, meta: dict, qc_dir: Path) -> dict:
    frame_ids = source.frame_indices(plate)
    positions = list(range(len(frame_ids)))
    if not positions:
        raise ValueError(f"No frames for plate {plate}")

    ref_pos = min(REFERENCE_FRAME_IDX, positions[-1])
    ref_full = source.read_frame(plate, ref_pos)
    cx, cy, r = detect_agar_circle(ref_full)

    # Which frame counts as the clean background?
    requested = meta.get("Baseline_Frame")
    requested = None if requested is None or pd.isna(requested) else int(requested)
    if requested is not None:
        baseline_pos = min(requested, positions[-1])
        baseline_source = "manual"
    elif AUTO_BASELINE:
        baseline_pos, _diag = find_settled_frame(source, plate, cx, cy, r)
        baseline_source = "auto"
    else:
        baseline_pos, baseline_source = 0, "default"

    offset_h = float(meta.get("Start_Offset_h") or 0.0)

    base_full = source.read_frame(plate, baseline_pos)
    base_crop, valid = crop_to_circle(base_full, cx, cy, r)

    # The detected circle is the dish rim, PLATE_DIAMETER_MM across, which fixes the
    # physical scale independently of how much of the disc we choose to analyse.
    mm_per_px = (PLATE_DIAMETER_MM / 2.0) / float(r)

    segmenter = PlateSegmenter(base_crop, valid)
    agar = segmenter.agar

    # Pass 1 intake: read each frame once, keep the colour metrics and the lightness
    # plane. Full frames are not retained — only what the second pass needs.
    used_positions = positions[baseline_pos:]

    # Drop an unusable tail if the workbook says so. A frame captured with the lamp
    # off or the tray being handled is black or blown out, and would otherwise be
    # taken as a huge change against the baseline.
    last_frame = meta.get("Last_Frame")
    if last_frame is not None and pd.notna(last_frame):
        kept = [p for p in used_positions if frame_ids[p] <= int(last_frame)]
        if len(kept) >= 2:
            dropped = len(used_positions) - len(kept)
            used_positions = kept
            if dropped:
                print(f"  {plate}: dropped {dropped} frame(s) after f{int(last_frame)}")
    colour_rows = []
    for pos in used_positions:
        frame = source.read_frame(plate, pos)
        crop, _ = crop_to_circle(frame, cx, cy, r)
        if crop.shape != base_crop.shape:
            crop = cv2.resize(crop, (base_crop.shape[1], base_crop.shape[0]))
        mr, mg, mb, lum = mean_colour(crop, agar)
        colour_rows.append((mr, mg, mb, lum))
        segmenter.add(crop)

    masks, seg_diags = segmenter.run()

    qc_wanted = {p for p in used_positions
                 if frame_ids[p] in QC_FRAME_INDICES} | {used_positions[0]}

    rows = []
    for k, pos in enumerate(used_positions):
        frame_idx = frame_ids[pos]
        colony = masks[k]
        seg_diag = seg_diags[k]
        mr, mg, mb, lum = colour_rows[k]
        am = area_metrics(colony, agar, mm_per_px)

        if SAVE_QC_FRAMES and pos in qc_wanted:
            qc_dir.mkdir(parents=True, exist_ok=True)
            label = f"{plate}_{meta.get('Role','?')}"
            if meta.get("Strain"):
                label += f"_{meta['Strain']}_r{meta.get('Replicate','')}"
            crop, _ = crop_to_circle(source.read_frame(plate, pos), cx, cy, r)
            if crop.shape != base_crop.shape:
                crop = cv2.resize(crop, (base_crop.shape[1], base_crop.shape[0]))
            cv2.imwrite(str(qc_dir / f"{label}_f{frame_idx:03d}.png"),
                        qc_overlay(crop, colony))

        rows.append({
            "Frame": frame_idx,
            # Time_h is measured from when the plates went into the RI92, so
            # experiments whose capture started late stay comparable.
            "Time_h": frame_idx * FRAME_INTERVAL_HOURS + offset_h,
            "Time_h_capture": frame_idx * FRAME_INTERVAL_HOURS,
            "R_mean": mr, "G_mean": mg, "B_mean": mb, "Luminance": lum,
            "Colony_pixels": am["colony_pixels"],
            "Area_percent": am["area_percent"],
            "Area_mm2": am["area_mm2"],
            "Radius_mm": am["equivalent_radius_mm"],
            "Radius_px": am["equivalent_radius_px"],
            "Medium_drift_L": seg_diag["drift"],
            "Seg_threshold": seg_diag["threshold"],
            "Medium_visible_frac": seg_diag["pool_frac"],
        })

    df = pd.DataFrame(rows).sort_values("Frame").reset_index(drop=True)
    df["Agar_center_x"] = cx
    df["Agar_center_y"] = cy
    df["Agar_radius_px"] = r
    df["Baseline_Frame_used"] = baseline_pos
    df["Baseline_source"] = baseline_source
    df["Start_Offset_h"] = offset_h
    df["Drift_frozen_at"] = segmenter.frozen_at if segmenter.frozen_at is not None else -1
    return {"table": df, "circle": (cx, cy, r), "baseline": baseline_pos,
            "frozen_at": segmenter.frozen_at}


# Workbook fields that change the computed numbers. Editing any of these must
# invalidate the cache; everything else is a label and can simply be refreshed.
CACHE_RELEVANT_FIELDS = ("Role", "Include", "Baseline_Frame", "Start_Offset_h", "Last_Frame")


def _settings_key(lookup: dict) -> str:
    """Short hash of the processing-relevant settings for one experiment."""
    parts = []
    for plate in sorted(lookup):
        meta = lookup[plate]
        parts.append(plate + "=" + "|".join(
            str(meta.get(f)) for f in CACHE_RELEVANT_FIELDS))
    return hashlib.md5(";".join(parts).encode()).hexdigest()[:8]


def process_experiment(source: PlateSource, experiment: str, lookup: dict) -> pd.DataFrame:
    # The cache key covers the workbook settings that change the numbers, not just
    # SEG_VERSION. Without this, editing Last_Frame (or a baseline, an offset, an
    # Include flag) silently has no effect: the run loads the table computed under
    # the old settings and never re-reads a single frame.
    cache_file = CACHE_DIR / f"{experiment}__seg{SEG_VERSION}_{_settings_key(lookup)}.pkl"
    qc_root = OUTPUT_DIR / "QC" / experiment

    # A cache hit skips process_plate entirely, so deleted QC overlays would never
    # come back. Reprocess when the overlays are wanted but missing.
    qc_missing = SAVE_QC_FRAMES and not any(qc_root.glob("*.png"))
    if USE_CACHE and cache_file.exists() and not qc_missing:
        print(f"[cache] {experiment}")
        with open(cache_file, "rb") as f:
            cached = pickle.load(f)
        # Refresh every label column from the current plate map, rather than only
        # filling in ones that are absent. A corrected strain name or role would
        # otherwise persist in the cache indefinitely. Settings that change the
        # numbers are handled by the cache key instead, so this is always safe.
        if not cached.empty:
            for col in PLATE_META_COLUMNS:
                cached[col] = cached["Plate"].map(
                    lambda p, c=col: lookup.get(p, {}).get(c))
        return cached
    if USE_CACHE and cache_file.exists() and qc_missing:
        print(f"[cache ignored] {experiment}: QC overlays missing, regenerating")
    frames = []

    for plate in source.plates():
        meta = lookup.get(plate)
        if meta is None:
            print(f"  {plate}: not in the plate map — skipped")
            continue
        if not meta.get("Include", True):
            print(f"  {plate}: Include=N — skipped")
            continue
        if meta.get("Role") in EXCLUDE_ROLES:
            print(f"  {plate}: {meta['Role']} — excluded by EXCLUDE_ROLES")
            continue

        result = process_plate(source, plate, meta, qc_root)
        df = result["table"]
        for key in PLATE_META_COLUMNS:
            df[key] = meta.get(key)
        df["Experiment"] = experiment
        frames.append(df)

        final = df.iloc[-1]
        base = result["baseline"]
        base_note = "" if base == 0 else f"  [baseline f{base} = {base * FRAME_INTERVAL_HOURS:.1f} h]"
        print(f"  {plate}: {meta['Role']:<22} {str(meta.get('Strain') or ''):<12}"
              f" final area {final['Area_mm2']:7.1f} mm²  r {final['Radius_mm']:5.1f} mm{base_note}")

    if not frames:
        return pd.DataFrame()

    out = pd.concat(frames, ignore_index=True)
    with open(cache_file, "wb") as f:
        pickle.dump(out, f)
    return out


# ================================================================ summaries
METRICS = ["Area_mm2", "Area_percent", "Radius_mm", "Luminance",
           "R_mean", "G_mean", "B_mean"]

# Plate-map columns copied onto every plate's time series. Keep this in one place:
# when a new grouping column was added to the workbook but not to this list, the
# whole run completed and then died in summarise() with KeyError: 'Condition'.
PLATE_META_COLUMNS = (
    "Plate", "Role", "Strain", "Replicate", "Spores",
    "Carbon", "Nitrogen", "Group", "Notes",
    "Condition", "Bio_Replicate", "Minus_Control_Omits", "Last_Frame",
    "Gelling_Agent",
)


def summarise(raw: pd.DataFrame):
    """
    Two summaries, because the replicates are nested.

    The 3 plates of a strain within one run are TECHNICAL replicates: one spore
    prep, one batch of plates, one machine run. Separate runs of the same condition
    are BIOLOGICAL replicates.

      run_mean  average of the technical replicates -> one curve per run
      nested    mean +- SD across runs, n = number of runs. The headline numbers:
                error bars reflect biological variation.
      pooled    every plate treated as independent, n = plates. Tighter error bars,
                but this is pseudoreplication and understates the true variance.
                Exported for comparison only.

    For a condition with a single run, nested and pooled differ only in n.
    """
    raw = raw.copy()
    # Runs share a 30 min grid; round so cross-run averaging lines up exactly.
    raw["Time_h"] = raw["Time_h"].round(3)

    # Restrict each group to the window all of its plates cover. Plates start at
    # different baseline frames and runs can end early (the agarose run stops at
    # 26 h), so outside the common window the average is taken over a changing set
    # of plates. That makes the mean step up or down as members enter and leave —
    # which shows up as physically impossible negative growth in the curve, even
    # though every individual plate is monotonic.
    span = raw.groupby(["Condition", "Group", "Experiment", "Plate"])["Time_h"].agg(["min", "max"])
    window = span.groupby(level=["Condition", "Group"]).agg(start=("min", "max"), end=("max", "min"))
    raw = raw.merge(window, left_on=["Condition", "Group"], right_index=True, how="left")
    dropped = ((raw["Time_h"] < raw["start"]) | (raw["Time_h"] > raw["end"])).sum()
    raw = raw[(raw["Time_h"] >= raw["start"]) & (raw["Time_h"] <= raw["end"])].drop(
        columns=["start", "end"])
    if dropped:
        print(f"  [summarise] trimmed {dropped} plate-frames outside each group's "
              f"common time window, so every group average has constant membership")

    keys = ["Condition", "Group", "Role", "Time_h"]

    run_mean = (raw.groupby(keys + ["Bio_Replicate", "Experiment"], dropna=False)[METRICS]
                .mean().reset_index())
    run_mean["n_plates"] = (raw.groupby(keys + ["Bio_Replicate", "Experiment"], dropna=False)
                            .size().values)

    g_nested = run_mean.groupby(keys, dropna=False)
    nested_mean = g_nested[METRICS].mean().reset_index()
    nested_std = g_nested[METRICS].std(ddof=1).reset_index()
    nested_mean["n_runs"] = g_nested.size().values
    nested_std["n_runs"] = g_nested.size().values

    g_pooled = raw.groupby(keys, dropna=False)
    pooled_mean = g_pooled[METRICS].mean().reset_index()
    pooled_std = g_pooled[METRICS].std(ddof=1).reset_index()
    pooled_mean["n_plates"] = g_pooled.size().values
    pooled_std["n_plates"] = g_pooled.size().values

    return {"run_mean": run_mean,
            "nested_mean": nested_mean, "nested_std": nested_std,
            "pooled_mean": pooled_mean, "pooled_std": pooled_std}


# Plot styling for the control roles. The two controls answer different questions and
# are easy to confuse at a glance, so they are given deliberately unlike marks rather
# than two similar dashed greys.
CONTROL_STYLE = {
    ROLE_CONTAM: dict(color="#8E44AD", linestyle="--", marker="^",
                      markersize=7, markevery=10, markeredgecolor="#4A235A"),
    ROLE_MINUS: dict(color="#111111", linestyle="-.", marker="s",
                     markersize=5, markevery=10),
    ROLE_WATER: dict(color="0.55", linestyle=":", marker=None),
}

X_AXIS_LABEL = "Incubation time (h)"


def condition_titles(raw: pd.DataFrame) -> dict:
    """
    Plot title for each condition: the variable actually under test — 'Urea',
    'Xylose', 'Ammonium sulfate' — taken from Minus_Control_Omits, which is exactly
    the component that condition varies.

    Where two conditions test the same component (ammonium sulfate on agar and on
    agarose) the gelling agent is appended, because otherwise two different figures
    would carry an identical title.
    """
    info = (raw.groupby("Condition")
            .agg(variable=("Minus_Control_Omits", "first"),
                 gel=("Gelling_Agent", "first")))
    titles, counts = {}, {}
    for cond, row in info.iterrows():
        name = (str(row["variable"]) or "").strip() or str(cond)
        counts[name] = counts.get(name, 0) + 1
    for cond, row in info.iterrows():
        name = (str(row["variable"]) or "").strip() or str(cond)
        name = name[0].upper() + name[1:] if name else str(cond)
        gel = (str(row["gel"]) or "").strip()
        if counts.get(name.lower() if False else name, 0) > 1 and gel:
            name = f"{name} ({gel})"
        titles[cond] = name
    return titles


def plot_condition(condition: str, summary: dict, metric: str, ylabel: str,
                   filename: str, title: str | None = None):
    """
    Mean curve per strain with an SD band, plus a faint line for each individual run
    so run-to-run variation stays visible and an outlier run is easy to spot.
    """
    m = summary["nested_mean"]
    s = summary["nested_std"]
    runs = summary["run_mean"]
    m = m[m["Condition"] == condition]
    s = s[s["Condition"] == condition]
    runs = runs[runs["Condition"] == condition]
    if m.empty:
        return

    n_runs = int(runs["Bio_Replicate"].nunique())

    plt.figure(figsize=(9, 6))

    # Draw samples first so the two controls sit on top and stay readable.
    order = sorted(m["Group"].unique(),
                   key=lambda g: m[m["Group"] == g]["Role"].iloc[0] != ROLE_SAMPLE)
    for group in order:
        gm = m[m["Group"] == group].sort_values("Time_h")
        gs = s[s["Group"] == group].sort_values("Time_h")
        role = gm["Role"].iloc[0]
        label = group.split(" | ", 1)[-1].replace("_", " ")
        style = dict(CONTROL_STYLE.get(role, {}))
        width = 2.6 if role in CONTROL_STYLE else 2.2

        line, = plt.plot(gm["Time_h"], gm[metric], label=label, linewidth=width,
                         zorder=4 if role in CONTROL_STYLE else 2, **style)

        if role == ROLE_SAMPLE and len(gs) == len(gm):
            sd = gs[metric].fillna(0).values
            plt.fill_between(gm["Time_h"], gm[metric].values - sd,
                             gm[metric].values + sd, alpha=0.16,
                             color=line.get_color(), linewidth=0)

        # Individual runs, faint.
        if n_runs > 1:
            for rep, gr in runs[runs["Group"] == group].groupby("Bio_Replicate"):
                gr = gr.sort_values("Time_h")
                plt.plot(gr["Time_h"], gr[metric], linewidth=0.9, alpha=0.45,
                         color=line.get_color(), linestyle="-")

    plt.xlabel(X_AXIS_LABEL)
    plt.ylabel(ylabel)
    # Title names only the variable under test, e.g. "Urea".
    plt.title(title or condition, fontsize=15)
    handles, labels = plt.gca().get_legend_handles_labels()
    if n_runs > 1:
        handles.append(plt.Line2D([], [], color="0.4", linewidth=0.9, alpha=0.6))
        labels.append("individual runs")
    plt.legend(handles, labels, fontsize=9)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    out = OUTPUT_DIR / "Plots" / f"{condition} — {filename}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=300)
    plt.close()


def growth_rates(raw: pd.DataFrame) -> pd.DataFrame:
    """Max radial extension rate (mm/h) per plate, from a rolling linear fit."""
    records = []
    for (exp, plate), g in raw.groupby(["Experiment", "Plate"]):
        g = g.sort_values("Time_h")
        t = g["Time_h"].values
        r = g["Radius_mm"].values
        window = max(3, int(6 / FRAME_INTERVAL_HOURS))   # ~6 h window
        best = 0.0
        best_t = np.nan
        for i in range(len(t) - window):
            slope = np.polyfit(t[i:i + window], r[i:i + window], 1)[0]
            if slope > best:
                best, best_t = slope, float(t[i])
        meta = g.iloc[0]
        records.append({
            "Condition": meta.get("Condition", exp),
            "Bio_Replicate": meta.get("Bio_Replicate", 1),
            "Experiment": exp, "Plate": plate, "Role": meta["Role"],
            "Strain": meta["Strain"], "Replicate": meta["Replicate"],
            "Group": meta["Group"],
            "Max_radial_rate_mm_per_h": best,
            "Window_start_h": best_t,
            "Baseline_Frame_used": int(meta["Baseline_Frame_used"]),
            "Baseline_source": meta["Baseline_source"],
            "Start_Offset_h": float(meta["Start_Offset_h"]),
            "Final_radius_mm": float(r[-1]),
            "Final_area_mm2": float(g["Area_mm2"].values[-1]),
        })
    return pd.DataFrame(records)


# ================================================================ main
def main():
    plate_map = load_plate_map(PLATEMAP_PATH)
    print(f"Plate map: {len(plate_map)} plates across "
          f"{plate_map['Experiment'].nunique()} experiment(s).\n")

    # Check now, not after an hour of image processing. Every column the summaries
    # group by has to exist on the plate map before any frame is read.
    absent = [c for c in PLATE_META_COLUMNS if c not in plate_map.columns]
    if absent:
        raise ValueError(
            f"Plate map is missing column(s) required by the analysis: {absent}. "
            f"PLATE_META_COLUMNS in this script and the workbook have diverged."
        )

    sources = discover_sources(INPUT_DIR)
    if not sources:
        raise ValueError(f"No plate data (folders or zips) found in {INPUT_DIR}")

    all_raw = []
    for source in sources:
        experiment = resolve_experiment(plate_map, source.name)
        if experiment is None:
            print(f"[SKIP] '{source.name}' has no row on the Experiments sheet "
                  f"of {PLATEMAP_PATH.name}")
            source.close()
            continue
        if ONLY_EXPERIMENTS and experiment not in ONLY_EXPERIMENTS:
            source.close()
            continue

        lookup = plate_lookup(plate_map, experiment)
        if not lookup:
            print(f"[SKIP] {experiment}: plate map sheet is empty")
            source.close()
            continue

        print(f"\n=== {experiment} ({len(source.plates())} plate folders) ===")
        raw = process_experiment(source, experiment, lookup)
        source.close()
        if not raw.empty:
            all_raw.append(raw)

    if not all_raw:
        raise ValueError(
            "Nothing was processed. Fill in the plate map sheets for the "
            "experiments you want to analyse, then rerun."
        )

    raw = pd.concat(all_raw, ignore_index=True)
    summary = summarise(raw)
    rates = growth_rates(raw)

    if EXPORT_CSV:
        raw.to_csv(OUTPUT_DIR / "raw_per_plate_timeseries.csv", index=False)
        # Headline: nested, error bars across biological replicates.
        summary["nested_mean"].to_csv(OUTPUT_DIR / "group_mean_timeseries.csv", index=False)
        summary["nested_std"].to_csv(OUTPUT_DIR / "group_std_timeseries.csv", index=False)
        summary["run_mean"].to_csv(OUTPUT_DIR / "run_mean_timeseries.csv", index=False)
        # For comparison only — see the docstring on summarise().
        summary["pooled_mean"].to_csv(OUTPUT_DIR / "pooled_mean_timeseries.csv", index=False)
        summary["pooled_std"].to_csv(OUTPUT_DIR / "pooled_std_timeseries.csv", index=False)
        rates.to_csv(OUTPUT_DIR / "growth_rates_per_plate.csv", index=False)

        # Growth rate per condition, both ways.
        samples = rates[rates["Role"] == ROLE_SAMPLE]
        if len(samples):
            per_run = (samples.groupby(["Condition", "Group", "Bio_Replicate"])
                       ["Max_radial_rate_mm_per_h"].mean().reset_index())
            nested_rate = (per_run.groupby(["Condition", "Group"])
                           ["Max_radial_rate_mm_per_h"]
                           .agg(mean="mean", sd="std", n_runs="count").reset_index())
            pooled_rate = (samples.groupby(["Condition", "Group"])
                           ["Max_radial_rate_mm_per_h"]
                           .agg(mean="mean", sd="std", n_plates="count").reset_index())
            nested_rate.to_csv(OUTPUT_DIR / "growth_rates_nested.csv", index=False)
            pooled_rate.to_csv(OUTPUT_DIR / "growth_rates_pooled.csv", index=False)

    with open(OUTPUT_DIR / "raw_data.pkl", "wb") as f:
        pickle.dump(raw, f)

    titles = condition_titles(raw)

    # ---- single-number metric and per-strain summary --------------------
    kr_rates = per_plate_rates(raw)
    minus_flags = minus_control_check(kr_rates)
    strain_summary = strain_substrate_summary(kr_rates, minus_flags)

    if EXPORT_CSV:
        kr_rates.to_csv(OUTPUT_DIR / "Kr_per_plate.csv", index=False)
        minus_flags.to_csv(OUTPUT_DIR / "minus_control_check.csv", index=False)
        strain_summary.to_csv(OUTPUT_DIR / "strain_substrate_summary.csv", index=False)
        if not strain_summary.empty:
            kr_matrix(strain_summary).to_csv(OUTPUT_DIR / "Kr_matrix.csv")
            # Exponential-phase equivalent, kept alongside the linear one.
            kr_matrix(strain_summary, value="Mu_per_h").to_csv(
                OUTPUT_DIR / "Mu_matrix.csv")

    if not strain_summary.empty:
        write_summary_workbook(strain_summary, kr_rates, minus_flags,
                               OUTPUT_DIR / "Strain_Summary.xlsx")
        if PLOTS:
            plots = OUTPUT_DIR / "Plots"
            # Two figures, split by the minus-control result, so the substrates you
            # can interpret are not read on the same colour scale as the ones you
            # cannot. Either is skipped if no substrate falls in that group.
            wrote_ok = plot_kr_heatmap(strain_summary, plots / "Kr heatmap.png",
                                       labels=titles, failed=False)
            wrote_bad = plot_kr_heatmap(
                strain_summary, plots / "Kr heatmap (minus control failed).png",
                labels=titles, failed=True)
            n_ok = int((~strain_summary["Minus_control_failed"].fillna(False)).sum())
            n_bad = int(strain_summary["Minus_control_failed"].fillna(False).sum())
            print(f"\nKr heatmaps: {n_ok} strain-substrate cells passed the minus-control "
                  f"check, {n_bad} failed")
            if not wrote_ok:
                print("  no substrate passed — 'Kr heatmap.png' not written")
            if not wrote_bad:
                print("  no substrate failed — 'Kr heatmap (minus control failed).png' "
                      "not written")

    if PLOTS:
        for condition in sorted(raw["Condition"].unique()):
            title = titles.get(condition, condition)
            plot_condition(condition, summary, "Radius_mm",
                           "Colony equivalent radius (mm)", "radial growth", title)
            plot_condition(condition, summary, "Area_mm2",
                           "Colony area (mm²)", "colony area", title)
            plot_condition(condition, summary, "Luminance",
                           "Mean luminance", "luminance", title)

    print("\nDone.")
    print(f"Outputs in: {OUTPUT_DIR}")

    pooled = raw.groupby("Condition")["Experiment"].nunique()
    multi = pooled[pooled > 1]
    if len(multi):
        print("\nPooled as biological replicates:")
        for cond, n in multi.items():
            runs = sorted(raw[raw["Condition"] == cond]["Experiment"].unique())
            print(f"  {cond}: {n} runs -> {', '.join(runs)}")

    if not strain_summary.empty:
        print("\nRadial extension rate Kr (mm/h), by strain and substrate:")
        m = kr_matrix(strain_summary)
        print(m.round(2).to_string())
        flagged = minus_flags[minus_flags["Minus_control_failed"]]
        if len(flagged):
            print("\nSubstrates to review — minus control grew nearly as well as the samples:")
            for _, row in flagged.sort_values("Minus_over_sample", ascending=False).iterrows():
                print(f"  {row['Condition'][:44]:46} minus/sample Kr = "
                      f"{row['Minus_over_sample']:.2f}  (omitted: {row['Omitted_component']})")
        print(f"\nStrain summary written to {OUTPUT_DIR / 'Strain_Summary.xlsx'}")


if __name__ == "__main__":
    main()
