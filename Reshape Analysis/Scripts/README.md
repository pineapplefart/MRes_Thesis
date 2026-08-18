# RESHAPE radial growth analysis

Pipeline for *Rhizopus microsporus* germling radial growth from RI92 timelapse
image sequences (30 min interval, 48 h, 15 plates per tray).

Adapted from Juan Sandoval's video-based RESHAPE pipeline in `../JuanCode/`.

## Files

| File | Purpose |
|---|---|
| `../RESHAPE_PlateMap.xlsx` | Which condition sits in each tray position, one sheet per experiment. **Fill this in first.** |
| `reshape_platemap.py` | Reads and validates the workbook. Run it alone to check your entries. |
| `reshape_growth_analysis.py` | The pipeline. Run this to produce curves, CSVs and QC images. |

## Tray layout

The RI92 tray is **3 rows x 5 columns**. In a plate code the letter is the column
(A–E) and the digit is the row (1–3):

|  | col A | col B | col C | col D | col E |
|---|---|---|---|---|---|
| **row 1** | A1 | B1 | C1 | D1 | E1 |
| **row 2** | A2 | B2 | C2 | D2 | E2 |
| **row 3** | A3 | B3 | C3 | D3 | E3 |

Each experiment sheet lists plates in that reading order, with a rule between rows,
so it matches a photo of the tray.

Nothing in the analysis depends on this — plates are keyed by the folder name the
RI92 wrote (`A1`, `B1`, …), and `Grid_Row`/`Grid_Col` exist only so you can sort and
filter by physical position.

## Workflow

1. Fill the `Strains` sheet in `RESHAPE_PlateMap.xlsx` with your 4 strain IDs.
2. For each experiment, fill the yellow columns from your phone photo of that
   tray. The sanity block under each table should read 15 / 12 / 1 / 1 / 1 / 12 / 12 / 0.
   Leave `Baseline_Frame` blank — see *Condensation* below.
3. Check the map parses cleanly:
   ```bash
   cd Scripts && python reshape_platemap.py
   ```
   Sheets you have not filled in yet are reported as "not filled in yet" and skipped,
   so you can start analysing before the whole workbook is done.
4. Run the analysis:
   ```bash
   python reshape_growth_analysis.py
   ```

## What changed from Juan's version

- **Image sequences, not video.** `PlateSource` walks `<Experiment>/<Plate>/*-PlateXX-N.jpg`
  and reads frames **straight out of the `.zip`** — no need to extract 23 GB.
  An extracted folder of the same name is used in preference if one exists.
- **No `x`/`y` centre columns.** RI92 exports each plate already cropped to its own
  3040x3040 frame, so the plate ID *is* the position (3 rows x 5 columns, see above). The agar circle is still fitted
  per plate with Hough, with a centred-circle fallback if detection fails.
- **Nothing intermediate is written to disk.** Frames are masked, cropped and
  segmented in memory. Only a handful of QC overlay PNGs per plate are saved.
- **`DOWNSCALE = 0.25`** (3040 -> 760 px). About 16x faster; area fractions shift by
  well under 1%. Set to `1.0` for full resolution.
- **Controls are first-class.** Water, contamination and minus-variable plates are
  labelled, plotted in grey dashed styles alongside the strains, and excluded from
  strain averaging.
- **`Include = N`** drops any plate (contaminated, fogged lid) without editing code.
- **Per-plate growth rates.** `growth_rates_per_plate.csv` gives the maximum radial
  extension rate (mm/h) from a rolling 6 h linear fit, plus final radius and area.
- **Automatic condensation handling** and **late-start correction**, below.

## Segmentation

Mycelium is brighter than the medium on every substrate here, so the pipeline
thresholds a signed rise in CIE lightness against the baseline frame rather than any
colour change. Three things make that hold up (`SEG_VERSION = 2`):

- **The dish wall is excluded** (`RIM_KEEP = 0.78`). The Hough circle lands on the
  outer rim; measured on water controls, the wall and meniscus artifact sits in the
  0.80–0.90 radial band. Analysing out to the rim caused the rim reflection itself to
  be selected as the colony. There is no agar beyond ~0.8, so nothing real is lost.
- **Medium drift is subtracted.** Opaque media darken markedly over 48 h — cellulose
  fell 22 L units. Each frame, the medium's own shift is measured from the pixels the
  colony has not yet claimed, so a global change cannot flood the plate.
- **Drift is frozen while the pool is still mostly medium** (`MIN_MEDIUM_FRAC = 0.60`).
  Once the colony dominates, "unclaimed" pixels are really colony; re-estimating from
  them walks the drift upward until the threshold exceeds the colony's own brightness
  and growth stops being detected. Molasses A1 drifted to +14 L and stalled at 58% on
  a fully covered plate. Freezing early fixes it — that plate now reads 89.8%.

  A cap on drift *rise* looks like the natural fix here and is the wrong one: the blue
  agars legitimately brighten, and capping put a 43% mask on a visibly clean Dextrose
  Urea plate at 24 h. The comment in the source says so; please don't reintroduce it.

### Two passes, not one — and why growth can never be negative

Segmentation runs **twice** over each plate (`PlateSegmenter.add` then `.run`):

1. Sweep forward, tracking medium drift and thresholding each frame. Take the frame
   where the selected colony is largest — normally near the end, where it is clearest
   — as the colony's **footprint**.
2. Replay every frame restricted to that footprint, accumulating. Because the mask is
   a union over time, it cannot shrink: **negative growth is impossible by
   construction**, not merely unlikely.

A single streaming pass cannot achieve this, and both obvious variants fail:

- *Re-pick the best blob each frame* → the choice jumps between objects and the mask
  shrinks. This produced the negative growth rates.
- *Lock onto the first blob and only grow it* → fixes shrinking, but freezes the plate
  on whatever it latched onto first. Molasses A1 locked onto an 8950 px artifact at
  frame 10 and stayed there for the remaining 43 h while the real colony, a separate
  component, was never picked up. Size cannot separate the two — that artifact was
  larger than many genuine early colonies.

**Exactly one region per plate.** All the spores go in as a single droplet, so there
is one colony. Each frame's mask is the connected component of the accumulated mask
that contains an anchor: the inoculum core on the first frames (the deepest part of
the footprint, by distance transform), then the previous frame's mask thereafter.

That single choice gives both guarantees at once. The previous mask is connected and
lies wholly inside the accumulated set, so the component holding it is a superset —
the colony can only grow, and it can never be reported as two pieces.

Getting there took four wrong turns, all recorded here so they are not retried:

| approach | outcome |
|---|---|
| union mask with the previous frame | monotone, but welds on a second disconnected region |
| anchor on a single centre point | sparse hyphae straddle it; plates read 0% while growth is visible |
| anchor on the first detection | that detection is sometimes a speck the colony never reaches, freezing the plate |
| anchor on the inoculum core throughout | keeps whichever component overlaps the anchor *most*, so a newly merged region can steal it and the mask shrinks |

Hole filling is also not quite monotone on its own — a hole can open onto the outside
as the colony grows — so the finished mask is fed back into the accumulator rather
than being recomputed from the raw detections each frame.

Two further details that matter more than they look:

- The footprint comes from **one frame**, not the union of all frames. Unioning 96
  frames of noise on a plate with no colony produced a large blob that everything
  then accumulated into — a blank cellulose control read 78% of the plate.
- During pass 1 only the *selected colony* is marked as claimed, never every
  thresholded pixel. Claiming raw detections lets scattered noise eat the medium pool
  on a blank plate, which freezes the drift early on a bad value and collapses the
  threshold.

**Choosing the colony.** Spores are pipetted into the middle, so candidate components
are scored on how close their *pixels* get to the plate centre (`NEAR_WEIGHT`,
`MEAN_WEIGHT`), and anything never reaching within `MAX_START_RADIUS` of the centre is
rejected outright. Scoring by centroid — the obvious approach, and what this did
originally — fails badly here: a long arc hugging the rim has a near-central centroid,
which is how rim reflections kept beating real colonies.

**One continuous colony per plate.** Loose aerial hyphae threshold as disconnected
filaments, so the mask is closed with a kernel of `GAP_CLOSE_FRAC` x plate radius to
bridge gaps between filaments of the same colony, then `FILL_HOLES` fills anything
enclosed by the outline. Closing is dilate-then-erode, so it fills gaps rather than
growing the colony outward. Selection is also re-run *after* the monotonic union —
without that, a component that won on an early frame survives forever as a detached
blob and the plate ends up with several "colonies".

This raises the measured area on ragged colonies, which is the intent: the fragmented
mask was measuring only the denser filaments inside the colony rather than its extent.
Dextrose Urea A1 at 48 h goes 12.6% -> 34.7%; already-solid colonies barely move.

`GAP_CLOSE_FRAC` is the sensitivity knob for loose aerial hyphae, and it is preferred
over lowering `MIN_DELTA_L`/`MAD_K` because it extends the mask along the colony
without lowering the bar for noise elsewhere. Measured at 48 h:

| `GAP_CLOSE_FRAC` | Urea A1 | Molasses A1 | blank plate |
|---|---|---|---|
| 0.035 | 31.0% | 86.7% | 0.6% |
| **0.06 (current)** | **34.7%** | **89.8%** | **0.7%** |
| 0.09 | 38.0% | 90.3% | 1.2% |

The faintest outer halo of aerial hyphae sits below threshold and is not captured at
any of these settings. Reaching it requires lowering `MIN_DELTA_L`/`MAD_K`, which
costs specificity badly — at `MIN_DELTA_L = 2.0` a visibly clean Urea plate read 43%
at 24 h and blanks grew phantom blobs up to 18%.

`Medium_drift_L`, `Seg_threshold` and `Medium_visible_frac` are written per frame to
`raw_per_plate_timeseries.csv` so any plate can be audited.

### Measured effect

| Plate | Before | After | |
|---|---|---|---|
| Cellulose D3 contamination ctrl | 5036 mm² | 501 mm² | flooding fixed; residual is a real colony (verified in the raw frames) |
| Cellulose C2 water | rim arc, 4.6% | 1.5% | rim artifact removed |
| Coffee C2 water | 553 mm² | 162 mm² | rim artifact removed |
| Molasses C2 water | 270 mm² | 83 mm² | rim artifact removed |
| Dextrose Urea A1 f71 | rim arc selected | colony selected | |
| Xylose A1 f71 sample | 3.4% | 14.4% | fluffy margin now caught |
| Molasses A1 f95 sample | flooded | tracks colony | |

Areas in mm² are calibrated from the detected dish rim (`PLATE_DIAMETER_MM` across),
not from the pixel count of the analysed region, so narrowing the disc did not change
the physical scale.

## Condensation

Growth is measured as the difference from a clean background frame. On plates that
start with condensation, frame 0 is not clean, and the condensation clearing gets
counted as a colony.

With `Baseline_Frame` left blank, the pipeline picks the background frame per plate.
It compares each of the first 8 frames against a settled reference at frame 20 (10 h —
past condensation, before growth), takes the plateau of that curve, and accepts the
first frame within 1.5x of it. A plate is only shifted off frame 0 if frame 0 is at
least 2x the plateau, so slow illumination drift doesn't trigger it. Frames before the
chosen baseline are dropped from that plate's series.

Measured on your data:

| Experiment | Detector |
|---|---|
| Cellulose AmmSulf | condensation on most plates, cleared by 1.0–2.5 h |
| Glycerol AmmSulf | milder condensation, cleared by 0.5–2.5 h |
| Coffee AmmSulf | clean, all plates keep frame 0 |
| Dextrose Urea | clean except C1 (2.0 h) |

Verified visually on Cellulose A1: a rim of droplets at frame 0, gone by frame 3,
which is what the detector chose.

The choice is recorded per plate as `Baseline_Frame_used` / `Baseline_source` in
`growth_rates_per_plate.csv`, and a QC overlay is always saved at the chosen frame.
To override, put a frame number in that plate's `Baseline_Frame` cell (1 frame = 30 min);
a value on the `Experiments` sheet sets a default for the whole experiment. The search
is capped at frame 8 (4 h), and growth doesn't start until ~14 h, so an over-cautious
baseline costs nothing.

## One number per plate: Kr

`reshape_summary.py` reduces each growth curve to the **radial extension rate Kr
(mm/h)** — the slope of colony radius against time, fitted over the linear phase only
(the frames where radius lies between 20% and 80% of that plate's final radius).

That window matters. It excludes the lag, which is a separate quantity that would drag
the slope down, and any plateau where the colony has run into the dish wall. So Kr
compares how *fast* a strain extends, independent of how long it took to start —
which is what you want when comparing one strain across substrates.

The same fit gives two things free:

- **Lag_h** — the fit's x-intercept, the extrapolated time at which growth began.
- **Kr_R2** — fit quality. On the current data Kr fits 129/131 sample plates with a
  median R² of 0.998. A few plates fit poorly (R² ~0.79) and those also produce
  nonsensical lags such as −47 h; **check R² before trusting a Kr**, since a curve
  that is not straight over the window has no single meaningful slope.

Kr is `NaN`, never 0, where a plate did not grow or the phase was too short to fit —
"no growth" and "not measurable" are different facts and stay distinct.

Averaging is nested, as elsewhere: technical replicates within a run first, then SD
across biological replicates.

### Exponential phase: µ

Alongside Kr, `Mu_per_h` is the **maximum specific growth rate** — the slope of
ln(radius) against time, taken from the steepest well-fitting 6 h window. Sliding the
window finds the exponential phase wherever it falls rather than assuming a position,
and the window's start time and R² are both reported so the choice is checkable.

Radii below `EXP_MIN_RADIUS_MM` (1 mm) are excluded from the fit. The logarithm
amplifies noise near zero, so without that floor the first speck to appear would
always look like the fastest growth on the plate.

**Kr and µ are not interchangeable.** Radial extension in filamentous fungi is
classically linear once the peripheral growth zone establishes; µ describes the
earlier expansion. On the current data they correlate at Pearson **r = 0.79** across
129 plates — related, but they rank plates differently often enough that the choice
matters. Both are kept for exactly that reason.

Measured: µ fits 129/131 sample plates, median R² 0.985, range 0.028–0.518 h⁻¹
(doubling times 1.3–24.5 h). The µ window starts at a median of 19 h, just after the
median Kr lag of 17.3 h — so it is picking up the phase immediately following the lag,
which is what it should be doing.

`Doubling_time_h` (ln2/µ) is included since it is often easier to reason about.

### Which substrates to exclude

A substrate whose **minus control grew nearly as well as its samples** did not
discriminate: the component under test made little difference, so its Kr says little
about that carbon or nitrogen source. Those are flagged wherever
`Minus_Kr / Sample_Kr >= 0.70` (`MINUS_FAIL_RATIO`), marked red in the workbook and
with `*` on the heatmap.

Outputs:

```
Strain_Summary.xlsx            README, Kr matrix, substrates to review, one sheet per strain
Kr_per_plate.csv               Kr, lag, R2, final radius for every plate
strain_substrate_summary.csv   Kr +- SD per strain per substrate, with flags
Kr_matrix.csv                  strain x substrate matrix, linear phase
Mu_matrix.csv                  strain x substrate matrix, exponential phase
minus_control_check.csv        minus-vs-sample ratio per condition, Kr and mu
Plots/Kr heatmap.png                        strain x substrate, minus control PASSED
Plots/Kr heatmap (minus control failed).png  strain x substrate, minus control FAILED
Plots/minus control comparison.png          minus controls against each other
Plots/linear Kr vs exponential mu.png       how far the two metrics agree
```

`Kr_per_plate.csv` carries both metrics per plate: `Kr_mm_per_h`, `Lag_h`, `Kr_R2`
for the linear phase, and `Mu_per_h`, `Doubling_time_h`, `Mu_window_start_h`, `Mu_R2`
for the exponential one.

### Comparing the minus controls

`Plots/minus control comparison.png` has two panels:

- **Left** — each minus control's Kr beside its own samples' median, so the controls
  can be ranked against one another in absolute terms.
- **Right** — each as a fraction of its own samples, with 1.0 and the flag threshold
  marked. This is the panel that answers whether omitting the component mattered.

Bars past 1.0 are the striking ones: the plate *without* the carbon or nitrogen source
outgrew the plates that had it.

The heatmap is split in two by the minus-control result. Keeping them apart matters
because each figure gets its **own colour scale** — on a combined figure the
uninterpretable substrates stretch the scale and flatten the contrast between the ones
you can actually use. Either figure is skipped, and nothing written, if no substrate
falls in that group.

## Adding a condition before the data exists

`Glutamate Ammonium Sulfate` is set up but empty, awaiting data. The pattern for any
future condition:

1. Add a row to the `Experiments` sheet (`Data_Folder`, `Sheet_Name`, sources,
   `Minus_Control_Omits`, `Gelling_Agent`).
2. Add a plate sheet named to match `Sheet_Name`.
3. Drop the zip into `Data/` when it exists.

A row whose data folder is not present is simply never reached, and an empty plate
sheet is reported as *"not filled in yet — skipped"* under notes rather than as an
error. So the workbook can run ahead of the data without breaking anything.

Two things to correct for glutamate once the run is done:

- **`Data_Folder` is a guess** — `Glutamate Ammonium Sulfate`, following the naming of
  the other runs. It must match the folder inside the zip exactly, so check it.
- **`Gelling_Agent` is blank and highlighted red**, as for the rest of the carbon
  series.

Worth noting for interpretation: glutamate is an amino acid, so it supplies **both
carbon and nitrogen**. Its sample plates carry glutamate *and* ammonium sulfate, while
the minus-glutamate control keeps the ammonium sulfate. The comparison is therefore
clean for carbon, but the samples have more available nitrogen than the control — so a
difference between them is not purely a carbon effect. The other carbon sources
(cellulose, xylose, fructose, glycerol, molasses, coffee) do not have this property.

## Incomplete trays and the minus control

- **`unused`** is a Role for a position deliberately not set up — a missing replicate,
  or a control that was not made. It is excluded from analysis and reported as a note
  rather than an error. Leaving the cell *blank* instead still warns, because blank is
  indistinguishable from "not filled in yet".
- Shortfalls against a full tray (12 sample / 1 water / 1 contamination / 1 minus) are
  printed under **"Plate map notes (not errors)"**. Currently: Cellulose is 11/12
  sample (no NH triplicate) and Dextrose AS V4 has no contamination control.
- **`Minus_Control_Omits`** on the `Experiments` sheet records what the minus control
  leaves out, and names it in plots — `Dextrose + Ammonium Sulfate | no ammonium
  sulfate` rather than the uninformative `minus_variable_control`. For every dextrose
  run the omitted component is the **nitrogen source**, not the dextrose; for the
  carbon series it is the carbon source.
- **QC overlays regenerate automatically** if deleted: a cache hit normally skips
  reprocessing, so the pipeline reprocesses an experiment whose QC folder is empty.

## Biological replicates

Runs sharing a `Condition` on the `Experiments` sheet are pooled. Blank means the run
stands alone. `Bio_Replicate` numbers the runs within a condition.

Currently:

| Data folder | Condition | Bio_Rep | Offset |
|---|---|---|---|
| Dextrose Ammonium Sulfate Fresh Spores | Dextrose + Ammonium Sulfate | 1 | 0 |
| Dextrose AS V4 | Dextrose + Ammonium Sulfate | 2 | 0 |
| Dextrose Ammonium Sulfate Agarose Part 2 | Dextrose + Ammonium Sulfate **(agarose)** | 1 | 9.8 |

### The agarose run — deliberately not pooled

The third dextrose/ammonium-sulfate run used **agarose**; the other two used agar.
Gelling agent is a real variable (agar carries sulphated galactans and trace nutrients
that agarose lacks), so pooling it would confound gelling agent with biological
variation. It is its own condition. Change the `Condition` cell to merge it if you
later decide the distinction does not matter.

It also lost power overnight and was restarted, giving 52 frames (26 h) rather than 96.
Its plates reached 2 mm equivalent radius at a median of 8.25 h versus 18.0 h on agar,
so they had roughly a 9.8 h head start before capture resumed — that is the
`Start_Offset_h`.

**That offset is estimated from a growth landmark, not measured.** It assumes agar and
agarose share the same intrinsic lag, which is exactly the kind of thing that might
differ between them. Consequences:

- Extension rate and final size from this run are usable.
- Lag time and any other time-anchored measure are **not** independent evidence — the
  alignment assumed them. Do not report a lag difference between agar and agarose
  from this run.
- Replace `Start_Offset_h` with the true outage duration if you can recover it from
  the RI92 logs, and the caveat disappears.

Also worth noting: the agarose plates were far more variable in reaching 2 mm
(4.0–14.5 h) than the agar plates (16–21 h).

**To add a third run:** drop the zip in `Data/`, add a row to `Experiments` (with the
same `Condition` string), and add a plate sheet named to match `Sheet_Name`. Positions
are randomised per run, so each run needs its own sheet.

### Nested vs pooled

The replicates are nested. The 3 plates of a strain within one run are **technical**
replicates — one spore prep, one batch of plates, one machine run. Separate runs are
**biological** replicates.

- **Nested (headline).** Average technical replicates within a run, then take mean ± SD
  across runs. `n` = number of runs. Error bars reflect biological variation.
  → `group_mean_timeseries.csv`, `group_std_timeseries.csv`, `growth_rates_nested.csv`
- **Pooled (comparison).** All plates treated as independent, `n` = plates. This is
  pseudoreplication: it counts technical replicates as independent and understates the
  true variance. Exported so you can compare, not for headline figures.
  → `pooled_mean_timeseries.csv`, `pooled_std_timeseries.csv`, `growth_rates_pooled.csv`
- **Per run.** One curve per run, used for the plot overlay. → `run_mean_timeseries.csv`

The difference is not cosmetic. On strain 6205 at 47.5 h across the two runs available:

| | value |
|---|---|
| mean | 3003 mm² (identical either way) |
| SD across runs (nested, n=2) | 370 mm² |
| SD across plates (pooled, n=6) | 539 mm² |
| run 1 mean / run 2 mean | 2741 / 3264 mm² |

The runs differ by ~19%, which is real biological variation that pooling would bury.

With only two runs the nested SD rests on two numbers and is weak — treat it as
indicative until the third run is in.

Plots show the mean per strain with an SD band and a faint line per run, so an outlier
run is visible rather than averaged away.

## Plot conventions

- **Title** names only the variable under test — `Urea`, `Xylose`, `Ammonium sulfate`
  — taken from `Minus_Control_Omits`, which is by definition the component that
  condition varies. Where two conditions test the same component (ammonium sulfate on
  agar and on agarose) the gelling agent is appended, so no two figures share a title.
- **X axis** is `Incubation time (h)`, measured from when the plates entered the RI92
  (see *Late starts*), not from the first captured frame.
- **The two controls are styled to be unmistakable**, since they answer different
  questions and were previously two similar dashed greys:

  | | style |
  |---|---|
  | contamination control | purple, dashed, triangle markers |
  | minus-variable control | black, dash-dot, square markers |
  | strains | solid coloured lines with an SD band |

  Controls are drawn on top of the sample curves so they stay readable, and the
  minus control is labelled by what it omits (`no urea`) rather than by its role.

`CONTROL_STYLE` and `X_AXIS_LABEL` near the plotting code are the single places to
change any of this.

## Trimming an unusable tail

`Last_Frame` on the `Experiments` sheet is the last frame index to analyse, inclusive.
Blank means use every frame. Set it where the end of a run is unusable — a lamp left
off or the tray handled during the final capture makes those frames black or blown
out.

Currently set for **Glutamate Ammonium Sulfate = 94**: frame 95 has a mean lightness
of 4 against ~136 on frame 94, and across plates it is either black or blown out white
with no dish visible.

Note the colony areas barely change when such a frame is dropped — a black frame
produces no positive lightness rise, so nothing is detected, and the accumulated mask
cannot shrink. What it protects is the **colour columns**: `Luminance`, `R_mean`,
`G_mean` and `B_mean` are plate averages taken directly from the frame, so a dead
frame lands in the CSVs and the luminance plots as a vertical drop.

## Late starts

`Start_Offset_h` on the `Experiments` sheet is the gap between the plates going into
the RI92 and the first captured frame. `Time_h` is measured from when the plates went
in, so experiments that started late still line up with the rest; `Time_h_capture`
keeps the raw frame time. Already filled in from your notes:

| Experiment | Offset |
|---|---|
| Dextrose Sodium Nitrate | 1.5 h |
| Coffee Ammonium Sulfate | 4.0 h |

Correct these or add others on the `Experiments` sheet. Plot titles note the offset
when one is applied.

## Outputs (written to `../Outputs/`)

```
raw_per_plate_timeseries.csv   every plate, every frame, all metrics
group_mean_timeseries.csv      NESTED mean by strain/control  <- headline
group_std_timeseries.csv       NESTED SD, across biological replicates
run_mean_timeseries.csv        one curve per run (technical replicates averaged)
pooled_mean_timeseries.csv     pooled mean, all plates as n  (comparison only)
pooled_std_timeseries.csv      pooled SD                     (comparison only)
growth_rates_per_plate.csv     max radial rate, final radius/area per plate
growth_rates_nested.csv        max radial rate per strain, +- SD across runs
growth_rates_pooled.csv        same, pooled across plates    (comparison only)
Plots/                         one set per CONDITION, pooled across its runs
QC/<experiment>/               overlay PNGs at frames 0, 24, 48, 71, 95
cache/<experiment>__segN.pkl   per-experiment results; N is SEG_VERSION
```

## Settings worth knowing

Near the top of `reshape_growth_analysis.py`:

- `ONLY_EXPERIMENTS` — restrict a run to named experiments. Empty means all.
- `EXCLUDE_ROLES` — roles skipped entirely: not processed, not plotted, not in any
  CSV. Currently `("water",)`. Set to `()` to bring water back. To drop a single
  plate instead, set `Include = N` on its row in the workbook.
- `USE_CACHE` — reuse `cache/*.pkl` from a previous run. Delete the pickle for an
  experiment after changing segmentation settings.
- `DIFF_THRESHOLD` (20) — colour distance from frame 0 that counts as mycelium.
  Raise it if condensation or lighting drift is being picked up as growth.
- `PLATE_DIAMETER_MM` (90) — the mm calibration. Hough locks onto the outer dish
  rim, so this should be the dish's outer diameter, not the agar diameter.
- `KEEP_ONLY_CENTRAL_OBJECT` (True) — keeps the single largest central component.
  Set False if you expect multiple colonies per plate.

## Known caveats

- Background subtraction assumes the plate does not move after the baseline frame.
  Condensation present at the start is handled; a plate knocked mid-run or a lid that
  fogs *later* is not. A step change or sudden drop in a curve usually means one of
  those happened — check the QC overlays for that plate before trusting it.
- `Radius_mm` is the **equivalent** radius, `sqrt(area/pi)`, not a measured leading
  edge. It is the right quantity for a roughly circular colony and understates the
  extent of an irregular one.
- The rolling-fit growth rate takes the maximum slope anywhere in the run, so on a
  plate with no real growth it reports noise rather than zero. Read it alongside
  `Final_area_mm2`.
