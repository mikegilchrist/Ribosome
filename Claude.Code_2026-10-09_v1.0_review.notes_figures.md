# Open review items: manuscript figures vs code

Written 2026-10-09 on obed, Claude Code session ef5348d2. Collected for later
review by MG and Ricardo; nothing here has been changed in the code or the
manuscript. Evidence and line numbers: `Claude.Code_2026-10-09_v1.0_code.map.md`
(Sections 2 and 5). Rmd = `Text/Paper final figures.Rmd`; tex =
`../polysomeModel/mRNA_Populations_PLOS_GENETICS.tex`.

## 1. Caption and code disagree (decide: regenerate or relabel)

1. **Fig 10, decapping rates.** Caption: mu = 2.2e-4, 1.7e-3, 5.7e-3. Code:
   panel A plots 2.22e-4, 3.30e-4, 2.31e-3 (rows 45001:50000, 40001:45000,
   20001:25000; Rmd L985-990); panel C uses 3.30e-4 (L1026-1029); panel D
   uses 2.31e-3 (L1045-1067). Panel B is correct. Figs 6-8 select the
   caption's values correctly (blocks 10, 6, 2), so the likely fix is to use
   those blocks. Changes a published figure: needs MG and Ricardo.
2. **Fig 11, statistic and gene count.** Caption: Spearman's rho, 850 genes,
   p < 1e-16. Code: `cor.test` with default method (Pearson) on log10 values
   (L1757), printed as rho (L1758); 704 rows (L1636, hard-coded `1:704`).
   Recompute with the data and correct whichever is wrong.
3. **Fig 7, imax.** Caption: imax = 41. Code: 40 (L132).
4. **Fig 5, mid initiation value.** Caption: 2e-2. Code: 0.01 (L242).

## 2. Provenance gaps (ask Ricardo)

5. **Fig 9 (`MRL.png`) has no source in the repo.** Closest code is the
   `TE_linear_*` plots (L1200-1281, block 10, imax 4/39/194), which cover the
   capped state only; the decapped and full panels have no source.
6. **Figs 2, 3, 4, 10 are assembled outside R** (`*_slices.png`,
   `parameter_histograms.png`, `Protein_Production_v2.png`); the assembly
   files and the Fig 2 histogram inputs (`../Data`, `../../Data`) are not in
   the repo.
7. **The Rmd cannot knit from a clean session.** Duplicate chunk label
   "Figure 4" (L1193, L1295); `Fsum_frac` used at L1106 but only assigned in
   commented-out code (L896); `text()` with no arguments at L1777. The
   published figures came from an interactive session.

## 3. Model questions

8. **Fig 10 normalisation sits at the grid edge.** Protein production rises
   monotonically in kappa, so "normalized to the maximal protein production"
   divides by the value at kappa = 0.5, the upper end of the chosen grid
   (Part A finding). MG recalls an analytic result bearing on this; check
   with Ricardo and in the manuscript.
9. **Top class under tau-only interference.** With `which.interfere =
   c("tau")` the code sets initiation to 0 at class imax and constant below
   it. MG: this is intended (2026-10-09). Confirm the manuscript states it.
10. **Duplicated yeast mu value.** Entries 3 and 4 of the Presnyak quantile
    list are identical (0.0038508177). Possibly a copy error for a missing
    quantile.

11. **Which dimensionless decapping rate to present.** The code already
    solves in dimensionless form (`rescale = T, rescale.val = "tau"`, used
    by every model figure). Open for the paper: express decapping as
    mu/tau, or as 9*mu*imax/tau (decapping rate times the time a ribosome
    takes to cross the transcript), which is the group that sets the
    decapped share of protein output (Part A identity). The choice changes
    how figures comparing imax values read. Handoff step C.6 (a separate
    dimensionless solver) was dropped on 2026-10-09 as redundant with the
    `rescale` option (MG).

## 4. Minor labelling (non-paper plots unless noted)

12. "log2 Scaled Initiation" axis on log10 data (Rmd L1253, 1264, 1276,
    1287); legend line types not matching plotted lines (L1284-1289,
    L1366-1369); quantile lines at 0.015/0.475/0.989 labelled 1/50/99%
    (`Parameter histograms.R` L132, 147). Relevant to the paper only if the
    `TE_linear_*` plots became Fig 9.
