# Code map: Ribosome package vs the polysomeModel manuscript

Written 2026-10-09 on obed, Claude Code session ef5348d2 (Part B of
`Claude.Code_2026-10-09_v1.0_handoff_test.first.refactor.md`). Code state:
commit 03be8d5 (solver and figure code unchanged since f70c00a). Manuscript:
`../polysomeModel/mRNA_Populations_PLOS_GENETICS.tex`. Nothing was deleted or
edited in Part B.

Method. Function definitions and the caller graph were extracted mechanically
(R `parse()` over every `R/*.R`; word-boundary grep for `name(` across all R
files and the Rmd, excluding the function's own body). Rmd chunk contents,
figure outputs, parameter blocks and the script files were read by haiku
subagents that returned tables; the claims that matter most (Fsum_frac,
Fig 11 correlation method, gene count, Fig 5 kappa values) were re-checked by
direct grep. "Dead" below means no caller outside its own file.

Conventions. Parameter matrix rows come from `CreateParMatrix` (expand.grid,
kappa fastest, 5000 kappa values `seq(1e-4, 0.5, 1e-4)`), so block k = rows
(k-1)*5000+1 : k*5000. The active yeast mu vector (Rmd L31, L149) is:

| block | rows | mu (1/s) | half-life (min) |
|---|---|---|---|
| 1 | 1:5000 | 0.01 | 1.2 ("sim") |
| 2 | 5001:10000 | 5.776e-3 | 2 |
| 3 | 10001:15000 | 3.851e-3 | 3 |
| 4 | 15001:20000 | 3.851e-3 | 3 (duplicate of block 3) |
| 5 | 20001:25000 | 2.310e-3 | 5 |
| 6 | 25001:30000 | 1.650e-3 | 7 |
| 7 | 30001:35000 | 8.252e-4 | 14 |
| 8 | 35001:40000 | 4.443e-4 | 26 |
| 9 | 40001:45000 | 3.301e-4 | 35 |
| 10 | 45001:50000 | 2.222e-4 | 52 |

Arabidopsis (`AT_` objects only, Rmd L131-145, imax 40): block 2 = 1.005e-3
(11.5 min), block 6 = 1.005e-4 (115 min), block 10 = 7.70e-6 (1500 min).

DII = density-independent initiation = `which.interfere = c("tau")`.
DDI = density-dependent initiation = `which.interfere = c("kappa","tau")`.
All model figures use lambda = 1, tau = 1, delta = 1e5.

## 1. Answers to the handoff's four questions

1. **Is anything in `Manuscript_Figures.R` not reproduced in the Rmd?**
   Yes, a lot, but none of it is in the paper. Of its 128 output filenames,
   15 are also written by the Rmd and 113 are not. The 113 are earlier-draft
   `Fig1_*`/`Fig2_*` exploratory plots (mean/SD/loading series over mu and
   imax), `contour_*`, `Variation_in_Ribosomal_loading_*`, and the four
   Arabidopsis histograms `Figure1_At_*` (L48-116). Of those, only the
   `Figure1_At_*` histograms are candidates for a paper panel (Fig 2C, 2D).
   Its 7 own functions (`cub1`, `cub2`, `pm_generate`, `calc_adjugate`,
   `calc_adjugate_w_mu`, `percent_var`, a copy of `filled.contour2`) have no
   callers outside the file.
2. **Is `Continuous approximation.R` referenced by the paper?** No. The paper
   says three times (tex L584, 627, 797) that the capped distribution "can be
   roughly approximated by a truncated Gaussian", but that is a verbal
   description. The file holds a continuum boundary-value problem (not
   Gaussian), discrete ODE chains, and one moment-matched `dtruncnorm` overlay
   used only for an interactive plot. It writes no files, does not source the
   solver, and cannot run top to bottom (uses `full.sys`, `trnorm` and
   `CalcUnmarkedClass_newtau` before they exist). A second, unrelated
   `dtruncnorm` kernel sits in `Parameter histograms.R` L345-347.
3. **What does `Reanalysis_20200914.R` do?** An Affymetrix microarray QC and
   normalisation reanalysis (RMA, GCRMA, MAS5 variants, MDS, heatmaps) of an
   Arabidopsis polysome-fraction dataset (TX/TL, fractions NP and 5-12, three
   replicates), plus an attempt to reconstruct TX from TL with dilution
   factors (10,000-permutation search, L975-1013). It never calls the solver,
   uses hard-coded `~/utk/Research/Project_1/...` paths, and depends on
   objects it never defines (`GDA`, `IND_*`, `col_m5calls_TX`, ...).
   Unrelated to this manuscript.
4. **Which Rmd chunk produced each manuscript figure?** Table in Section 2.
   Fig 10 is mislabeled as the review found, and the mislabel is wider than
   the review said: three of its four panels use the wrong mu. Figs 3-8 use the
   mu values their captions state. Figs 5, 7 and 11 have other caption/code
   disagreements, listed in Section 5.

## 2. Manuscript figure -> code

| Fig | Paper image | Producing code | Blocks used (mu) | Caption says | Agrees? |
|---|---|---|---|---|---|
| 1 | `Figure1_biomodel_V3.png` | none (schematic) | - | - | n/a |
| 2 | `2023-07-04_parameter_histograms.png` | not in repo; panels likely from `Parameter histograms.R` (yeast length L72, half-life L131), `Manuscript_Figures.R` L48-116 (Arabidopsis length, half-life), `scaled initiation histogram.R` L20 (panel E), `Parameter histograms.R` L294/L326 (panel F, log10 FC); assembly not in repo | data | data | cannot check without `../Data` |
| 3 | `2023-07-04_Unmarked_slices.png` | panels: Rmd "Figure 5" chunk L1496 (`_Unmarked_DDI_Distributions_<i>`), heatmap Rmd "Figure 3" chunk L785 (`Figure1_Unmarked_pct`); assembly not in repo | block 10 (2.22e-4); slices at local rows 2, 30, 1030, 2030 = kappa 2e-4, 3e-3, 0.103, 0.203 | mu 2e-4; kappa 2e-1, 1.03e-1, 3e-3, 2e-4 | yes |
| 4 | `2023-07-04_Marked_slices.png` | panels: Rmd L1515 (`_Marked_DDI_Distributions_<i>`), heatmap L809 (`Figure1_Marked_pct`); assembly not in repo | block 10 | as Fig 3 | yes |
| 5 | `2023-07-09_Figure1_DIIvsDDI_medianlength_low_marking_with_labels.png` | Rmd "Figure1 1" chunk L241 (labels added outside R) | block 10; kappa `c(0.001, 0.01, 0.1)` (L242) | mu 2e-4; kappa 1e-3, **2e-2**, 1e-1 | mu yes; mid kappa no (0.01 vs 0.02) |
| 6 | `2023-07-09_Figure2_Marking_Rate_range_medianlength_with_labels.png` | Rmd "Figure 2" chunk L401 | blocks 10, 6, 2 (2.22e-4, 1.65e-3, 5.78e-3) | 2e-4, 1.7e-3, 5.7e-3 | yes |
| 7 | `2023-07-09_Figure2_At_Marking_Rate_range_medianlength_with_labels.png` | Rmd "Figure 2" chunk L554 | AT blocks 10, 6, 2 (7.7e-6, 1.0e-4, 1.0e-3); **imax 40** (L132) | 7.7e-6, 1e-4, 1e-3; **imax = 41** | mu yes; imax no |
| 8 | `2023-07-28_logodds.png` | Rmd "Figure 5" chunk L1432 | yeast and AT blocks 10, 6, 2 | yeast 2e-4, 1.7e-3, 5.7e-3; AT 7.7e-6, 1e-4, 1e-3 | yes |
| 9 | `MRL.png` | **not found**: no `MRL` or matching filename anywhere. Closest: Rmd "Figure 4" chunk L1200-1281 (`TE_linear_*`), mean ribosome count vs kappa at imax 4/39/194, block 10, but capped state only | block 10 | imax 4, 39, 194; mu 2.2e-4; panels capped, decapped, full | panel A plausibly; B and C have no source |
| 10 | `2023-07-17_Protein_Production_v2.png` | Rmd "Figure 3" chunk: panel A L985 (`protein_production_stacked_log_correct`), B L1006 (`_low_marking`), C L1026 (`_mid_marking`), D L1045 (`_high_marking`); imax 39; `_v2` assembly not in repo | A: blocks 10, **9**, **5**; B: 10; C: **9** (3.30e-4); D: **5** (2.31e-3) | 2.2e-4, 1.7e-3, 5.7e-3 | **no**: A (2 of 3 curves), C, D use the wrong blocks; under the convention of Figs 6-8 they should be blocks 6 and 2 |
| 11 | `Duc_Song_vs_model_log.png` | Rmd "Duc and Song" chunk L1751 | per-gene parameters, not blocks; **704 genes** (L1636) | **850 genes**; **Spearman's rho**, p < 1e-16 | **no**: code uses `cor.test` default (Pearson, L1757) on log10 values, labelled rho; gene count differs |

Fig 10 note. Figs 6-8 select blocks by value (`mu == par.ranges[[2]][k]`,
k = 10, 6, 2); Fig 10 alone uses hard-coded row ranges (L986-990,
L1007-1067). If the duplicated 3rd/4th mu entry were removed, rows
20001:25000 would hold 1.65e-3 (the caption's "median"), so the hard-coded
ranges may predate the duplicate. That is a hypothesis only: no nine-value
yeast mu list survives in the repo (`Manuscript_Figures.R` uses only the
Arabidopsis list), and it would not explain panel D's "high" label. Ask
Ricardo.

## 3. Functions

Status: live = reached from the Rmd; superseded = has a live replacement;
dead = no caller outside its own file (or none at all); unknown = needs MG.

### 3.1 `R/Polysome_functions_RAUC.R` (685 lines; the only file the Rmd sources, L10)

| Function | Lines | Status | Computes | Feeds | Manuscript equation |
|---|---|---|---|---|---|
| `CreateParMatrix` | 7-36 | live | expand.grid parameter matrix, kappa fastest | all model figs | - |
| `CalcUnmarkedClass` | 41-80 | live | capped steady state, tridiagonal solve (limSolve) | all model figs | eq:kappa_i, eq:Capped_ODE, eq:capped_solution |
| `CalcMarkedClass` | 82-131 | live | decapped steady state by tail sums | all model figs | eq:Decapped_ODE, eq:S, eq:decapped_solution |
| `Classes` | 134-153 | live | capped + decapped | Figs 5-11 | eq:decapped_distribution |
| `CHECK_STEADY_STATE` | 158-228 | dead, buggy | residual check; uses tau_i = i and returns after class 0 (KNOWN_FAILURES #4, #5) | none | - |
| `plot.series`, `plot.bar` | 237-267 | dead | base plots | none (callers only in samples.R) | - |
| `plot.bars`, `plot.bars.mat` | 270-371 | dead | bar plots | none (samples.R, Model Figures.R) | - |
| `series.mat.plot` | 373-393 | dead | PNG per parameter series | none (Model Figures.R) | - |
| `s2_approx` | 397-406 | dead | analytic variance approximation | none (Model Figures.R) | - |
| `CalcUnmarkedClass_newtau` | 412-450 | superseded | identical to base + `print()` | none (Continuous approximation.R only) | - |
| `CalcMarkedClass_newtau` | 455-504 | superseded | identical; falls back to base `CalcUnmarkedClass` | none | - |
| `Classes_newtau` | 507-526 | superseded | identical | none | - |
| `Ribo_percent` | 531-538 | live | column-normalised probabilities | Figs 3-8 | - |
| `Ribo_mean` | 540-548 | live | mean ribosomes per mRNA | Figs 9-11 | eq:Expected_ribo_load, eq:MRL |
| `Ribo_sd` | 551-565 | live (matrix); vector branch broken | SD of class | DII/DDI summaries | - |
| `Ribo_mode` | 567-576 | live | modal class | DII/DDI summaries | - |
| `Ribo_frac` | 579-581 | live | ratio | `M_ribo_frac_*` (not in paper) | - |
| `Ribo_frac_load` | 584-587 | live | load-weighted ratio | `M_ribo_frac_load_*` (not in paper) | - |
| `Check_dir` | 593-608 | live | create/set output dir | Rmd, Parameter histograms.R | - |
| `filled.contour2` | 612-684 | live | contour plot helper | Rmd "Figure1 1" | - |

Protein production (Fig 10) is not a function: Rmd L977-983 compute
`colSums(classes) * Ribo_mean(classes)` normalised to its max (eq:System_ribo_load).

### 3.2 Other files

| Function | File | Lines | Status | Computes |
|---|---|---|---|---|
| `CalcUnmarkedClassNP`, `CalcMarkedClassNP`, `ClassesNP` | samples.R | 4-47 | superseded | base solver with tau_NP = tau/(9 imax), always kappa+tau interference (equivalence proven in tests, KNOWN_FAILURES #2) |
| `plot.series`, `plot.bar`, `plot.bars`, `plot.bars.mat` | samples.R | 83-219 | dead | duplicates of RAUC plot helpers |
| `MakeUnmarkedMatrix`, `SIM_UMAT`, `MakeMarkedMatrix`, `SIM_MMAT` | samples.R | 474-640 | dead | generator matrices, 1.5M-step Euler simulation |
| `f` | samples.R | 456 | dead | one-line helper |
| `cont`, `marked.cont`, `unmarked.full`, `marked.full`, `full.sys`, `s2` | Continuous approximation.R | 17-310 | dead | continuum BVP, ODE chains, variance expression |
| `Ribo_percent`, `Ribo_mean`, `series.mat.plot` | Model Figures.R | 828, 1001, 1189 | dead | local redefinitions of RAUC functions |
| `CalcUnmarkedAdj_noM` | Model Figures.R | 871-918 | dead | capped steady state via adjugate |
| `cub1`, `cub2`, `pm_generate`, `calc_adjugate`, `calc_adjugate_w_mu`, `percent_var`, `filled.contour2` | Manuscript_Figures.R | 2425-3319 | dead | exploratory helpers; adjugate solvers |
| `harmonic` | Parameter histograms.R | 360-366 | dead | toy expected-load calc |
| `affy.2.tair`, `cplot`, `cplot2` | Reanalysis_20200914.R | 47-122 | dead | probe annotation, scatter plots |
| `hello` | hello.R | 16-18 | dead | template |

Called but defined nowhere in `R/`: `CalcUnmarkedClassNP_I`, `_NI`
(samples.R L68, 73), `CalcUnmarkedClassNPM`, `CalcMarkedClassNPM`
(Model Figures.R L231-264).

## 4. Rmd chunks (`Text/Paper final figures.Rmd`)

| Chunk | Lines | Status | Computes | Feeds |
|---|---|---|---|---|
| setup | 9-14 | live | sources RAUC; viridis, GenomicFeatures | all |
| DII | 19-72 | live | DII classes and summaries, imax 39 | Fig 5 |
| DII imax | 76-118 | live | DII over imax 4/12/39/98/194 | Fig 9 candidate, `M_ribo_frac_*_DII` |
| DDI | 123-181 | live | AT DDI at imax 40; yeast DDI at imax 39 | Figs 3-8, 10 |
| DDI imax | 185-221 | live | DDI over imax range | Fig 9 candidate, protein production at imax 194 |
| create directory | 224-229 | live | `setwd("../../Results/<date>")` | output location |
| Figure1 1 | 235-390 | live | DII vs DDI profiles and heatmaps | **Fig 5** |
| Figure 2 | 395-709 | live | marking-rate heatmaps, yeast and AT | **Figs 6, 7** |
| Figure 3 | 715-1185 | live, **does not run clean** | heatmaps, protein production, colour keys | **Figs 3B, 4B, 10** |
| Figure 4 | 1193-1292 | live | `TE_linear_*`, mean load vs kappa by imax | Fig 9 candidate |
| Figure 4 (duplicate label) | 1295-1316 | unknown | relative mRNA vs mu (no model output) | none in paper |
| Figure 5 | 1319-1566 | live | log-odds; per-kappa distribution panels | **Figs 3A, 4A, 8** |
| Duc and Song model results | 1569-1791 | live, needs `../../Data` | per-gene prediction vs Weinberg 2016 | **Fig 11** |

The duplicate chunk label "Figure 4" stops `knitr` from knitting the file at
all ("duplicate label"), and chunk labels do not match paper figure numbers.

## 5. Findings for MG and Ricardo (no code changed)

1. **Fig 10 caption vs code** (Section 2). Panels A, C, D plot mu 3.3e-4 and
   2.3e-3 under labels 1.7e-3 and 5.7e-3. Legends in code (L937, L971, L999)
   carry the same wrong labels. Needs a decision: regenerate with blocks 6
   and 2, or relabel.
2. **Fig 11 statistic.** Caption: Spearman's rho on 850 genes. Code: Pearson
   r on log10 values (L1757), labelled rho, 704 rows (L1636 hard-codes
   `1:704`). The paper's "rho" and gene count should be checked against data.
3. **Fig 7 imax.** Caption 41, code 40.
4. **Fig 5 mid kappa.** Caption 2e-2, code 1e-2.
5. **Fig 9 has no source in the repo**, and Figs 2, 3, 4, 10 are assembled
   (or labelled) outside R. Ricardo likely has the assembly files.
6. **The Rmd cannot knit from scratch.** Duplicate chunk label (Section 4);
   `Fsum_frac` used at L1106 but only assigned in commented-out code (L896),
   so the "Figure 3" chunk errors at that line in a clean session;
   L1777 `text()` with no arguments errors before `dev.off()`. The published
   figures were therefore made from an interactive session whose state is not
   fully reproducible from this file.
7. **Fig 10 normalisation** ("normalized to the maximal protein production
   across all parameters") divides by a maximum that sits at the kappa = 0.5
   edge of the grid (Part A finding). MG recalls a relevant analytic result.
8. **Labels**: "log2 Scaled Initiation" axis on log10 data (L1253, 1264,
   1276, 1287); legend line types not matching the plotted lines (L1284-1289,
   L1366-1369); quantile lines at 0.015/0.475/0.989 labelled 1/50/99%
   (`Parameter histograms.R` L132, 147). None of these are paper figures
   except possibly Fig 9.

## 6. What this means for Part C

- **The Rmd sources only `Polysome_functions_RAUC.R`.** Every other file in
  `R/` can leave `R/` without touching the figure source. Suggested home:
  `inst/scripts/` (or a top-level `Analyses/`), keeping them in git.
- **Keep in `R/`:** the live rows of 3.1. **Remove from RAUC:** `_newtau`
  set, plot helpers, `series.mat.plot`, `s2_approx` (move to scripts if
  wanted), and `CHECK_STEADY_STATE` (broken; the test suite replaces it).
- **`samples.R` NP set:** redundant (test-proven mapping); move with
  `samples.R` to scripts, do not port.
- **Dead files:** `hello.R`, `polysome_model_results.Rmd` (templates).
- **`Reanalysis_20200914.R`:** unrelated analysis; move out (or to its own
  repo), MG's call.
- Once `R/` holds only function definitions and `DESCRIPTION` is fixed,
  `R CMD check` becomes meaningful; that is the natural first Part C commit
  after the template removal.
- Fixing the Rmd (findings 1-4, 6) changes paper figures and needs MG's and
  Ricardo's sign-off; it is not a refactor.
