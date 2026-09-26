# GSP in Oncology · Data workbench

An open Shiny application for the **Good Statistical Practice in Oncology** hands-on workshop.
Participants upload a trial dataset and reproduce every analysis taught in the sessions,
from data checks to the result sentence.

## What it does

| Tab | Analyses |
|---|---|
| 1 · Data | Upload (CSV, TSV, Excel), variable mapping, six data checks, Table 1 by arm |
| 2 · Survival | Kaplan–Meier curves with numbers at risk and CI bands, medians, landmark rates, median follow-up by reverse Kaplan–Meier, Kaplan–Meier by hand |
| 3 · Compare arms | Groups by arm or any baseline factor (alone or combined), log-rank, Cox model with optional adjustment, proportional-hazards check, RMST, landmark difference and NNT |
| 4 · Subgroups | Forest plot of hazard ratios with interaction tests |
| 5 · Sample size | Time-to-event (events and patients), non-inferiority, two proportions, two means, Simon two-stage and A'Hern |
| 6 · Diagnostics | 2 × 2 table, sensitivity, specificity, predictive values, likelihood ratios, predictive values across prevalence, ROC with AUC and Youden cut-off |
| 7 · Agreement | Bland–Altman, ICC(2,1) and ICC(3,1), Cohen's kappa |
| 8 · Report | Result sentences and a CSV of the key results |

Display settings on the Data tab convert time between days, weeks, months and years, and
set group colours from a palette catalogue or a custom colour picker.

## Try it

Click **Use the example dataset** in the app, or upload `data/gsp_exercise_dataset.csv`
(or the `.xlsx`, which also carries a data dictionary). Both are simulated; no patient data
are held in this repository or by the app.

## Run locally

```r
install.packages(c("shiny", "bslib", "DT", "readxl", "survival", "ggplot2"))
shiny::runApp("app.R")
```

## Deploy to Posit Connect Cloud

1. Keep `app.R` and `manifest.json` together at the top of the repository. The manifest
   records the R version, the package versions and a checksum of `app.R`.
2. At connect.posit.cloud choose **Publish → Shiny**, select this repository and branch,
   and confirm `app.R` as the primary file.
3. Open **Share** on the published content and set **Public access** to *Enabled*.
4. Whenever `app.R` or the package set changes, regenerate the manifest with
   `write_manifest.R` (kept on your computer, outside the repository) and commit both files together.

## Licence

MIT — see `LICENSE`.
