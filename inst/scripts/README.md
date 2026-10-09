# inst/scripts: analysis scripts, not package code

Moved here from `R/` on 2026-10-09 (Part C.2). Building a package runs every
file in `R/` top to bottom, and these files are scripts: they read data from
`../Data` or `../../Data`, call `setwd()`, and write figures. With them in
`R/` the package could not be installed or checked. Contents are unchanged;
see `Claude.Code_2026-10-09_v1.0_code.map.md` at the repo root for what each
one does.

The manuscript figures come from `Text/Paper final figures.Rmd`, which uses
only `R/Polysome_functions_RAUC.R`. None of these scripts produce a
manuscript figure, with two possible exceptions noted in the map (Fig 2
panels from `Parameter histograms.R`, `Manuscript_Figures.R` and
`scaled initiation histogram.R`).

| File | What it is |
|---|---|
| `Manuscript_Figures.R` | Earlier version of the figure Rmd; mostly superseded |
| `Model Figures.R` | Older exploratory model figures |
| `Parameter histograms.R` | Yeast parameter histograms from data (Fig 2 candidates) |
| `scaled initiation histogram.R` | Duc et al. 2018 scaled initiation histogram (Fig 2E candidate) |
| `Continuous approximation.R` | Continuum BVP and ODE checks; not used by the paper |
| `samples.R` | `*NP` solver variants (redundant with the main solver; see tests) and simulations |
| `Reanalysis_20200914.R` | Arabidopsis microarray QC; unrelated to this manuscript |

Paths inside the scripts are relative to the working directory they were
run from, not to this folder, so moving them does not change how they run.
