# Handoff: Annotated notebooks + gene-count graphs

**Date:** 2026-08-17
**Task:** Analyze the existing genome-annotation-enhancement and promoter-panel code, and produce graphs showing how many genes/proteins survived each step (analyzed → matched → structure found → characterized/uncharacterized, etc.). Original notebooks were left untouched; annotated copies were created instead.

## What was produced

Three new notebooks, each a copy of an existing one with markdown explanations added above every step and executed graph cells appended at the end. No analysis logic was changed except one bug fix (noted below).

| New file | Location | Copy of | Status |
|---|---|---|---|
| `genomeAlignment_annotated.ipynb` | `mmseq2 alignments/` | `genomeAlignment.ipynb` | Executed, 4 charts |
| `newsamplingmethod_annotated.ipynb` | `Promoter_panel/` | `newsamplingmethod.ipynb` | Executed, 3 charts |
| `findingPromoters_annotated.ipynb` | `Promoter_panel/` | `findingPromoters.ipynb` | Executed, 1 chart |

All three run end-to-end in Jupyter with no errors. Chart PNGs were also saved alongside each notebook (same folder) so they can be pulled into slides/reports without reopening Jupyter.

## Why three notebooks, and how they connect

There are two separate pipelines in your code, and the second depends on the first:

1. **`genomeAlignment.ipynb`** (`mmseq2 alignments/`) — cross-annotates your BRAKER-predicted gene isoforms against four reference proteomes (*R. oligosporus*, ATCC 11559, ATCC 52813, *R. delemar*) via MMseqs2, filters to confident matches, checks AlphaFold for a predicted structure, and classifies each match as Characterized / Uncharacterized / Not Found. It exports `52813_dictionary.csv`, a list of isoforms with a confirmed structure and a real (non-"Uncharacterized") description.
2. **`findingPromoters.ipynb`** and **`newsamplingmethod.ipynb`** (`Promoter_panel/`) — pick a 12-gene promoter panel from your RNA-seq data. `newsamplingmethod.ipynb` is the **current** method and is the one that actually reads `52813_dictionary.csv` from step 1 to restrict candidates to functionally-annotated genes. `findingPromoters.ipynb` is an **earlier, superseded** method (simpler median-expression + annotation-count filter) kept for the record — it does not use the AlphaFold dictionary.

## Key numbers (from the executed notebooks)

**Genome annotation (`genomeAlignment_annotated.ipynb`)**
- Total BRAKER-predicted isoforms: **27,218**
- Matched to reference proteome at ≥95% identity: ATCC 52813 = 6,504; ATCC 11559 = 9,224; R. oligosporus and R. delemar also computed
- Of matches checked against AlphaFold — ATCC 52813: 8,202 Characterized / 1,003 Uncharacterized / 195 Not Found; ATCC 11559: 9,250 Characterized / 1,766 Uncharacterized / 214 Not Found
- Final coverage: **10,648 of 27,218** isoforms (39%) ended up with some cross-strain annotation; 16,570 remain unannotated

**Promoter panel, current method (`newsamplingmethod_annotated.ipynb`)**
- Genes with RNA-seq data: **21,522**
- AlphaFold-characterized genes (from step 1's dictionary): **7,517**
- Overlap (characterized genes with expression data — the real candidate pool): **7,425**
- Final panel: **12 genes**, chosen to be evenly spaced across the expression×stability Score distribution of that candidate pool

**Promoter panel, earlier method (`findingPromoters_annotated.ipynb`)**
- Genes with RNA-seq data: 21,522 → genes with ≥5/6 BRAKER annotations: 6,634 → overlap: 5,771 → final panel: 12 genes

## Note on the one code change

`findingPromoters.ipynb`'s last cell called `pd.to_csv(promoter_set, "promoter_set.csv")`, which errors (should be a method call on the DataFrame, not on the `pandas` module). This crashed execution before the graph cell could run, so in the **annotated copy only** it was fixed to `promoter_set.to_csv("promoter_set.csv")`. The original file was not touched — you may want to apply the same one-line fix there if you still use that notebook.

## Where things are

- `mmseq2 alignments/genomeAlignment_annotated.ipynb` + 4 PNGs (`genes_matched_per_strain.png`, `funnel_total_to_alphafold_checked.png`, `classification_characterized_vs_uncharacterized.png`, `final_annotation_coverage.png`)
- `Promoter_panel/newsamplingmethod_annotated.ipynb` + 2 PNGs (`promoter_selection_funnel.png`, `final_panel_score_distribution.png`)
- `Promoter_panel/findingPromoters_annotated.ipynb` + 1 PNG (`promoter_selection_funnel_v1.png`)

## Suggested next steps

- If you present this work, the funnel chart in `genomeAlignment_annotated.ipynb` (Step 2 cell) and the funnel in `newsamplingmethod_annotated.ipynb` together tell the full story: predicted isoforms → cross-strain annotation → promoter candidate pool → final panel.
- Consider deciding whether `findingPromoters.ipynb` should be archived/deleted now that `newsamplingmethod.ipynb` is the adopted method — it's superseded but still runnable.
