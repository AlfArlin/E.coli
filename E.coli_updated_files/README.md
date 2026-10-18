# Proteome allocation and growth rate in *E. coli* K-12

**Author:** Arlin Alfaro Arias · **Course:** MSB1015 Scientific Programming, Maastricht University · October 2026

## Research question

Which functional protein groups (COG categories) of the *E. coli* K-12 proteome take a larger or
smaller share of the proteome as growth rate changes across nutrient conditions, and does that
allocation pattern hold across independent studies?

## Data

Four published, **already-processed** absolute proteomics datasets are combined. No raw mass
spectrometry data are reprocessed: each study is used in the processed form its authors published,
and the pipeline only harmonises units, identifiers, annotation and metadata.

| label | study | technology / original unit | files in `data/raw/` |
|---|---|---|---|
| DF01 | Schmidt A. *et al.* (2016) The quantitative and condition-dependent *Escherichia coli* proteome. *Nat. Biotechnol.* 34:104–110. [doi:10.1038/nbt.3418](https://doi.org/10.1038/nbt.3418) | DDA, protein copies per cell | `1.1The quantitative and condition-dependent Escherichia coli proteome.xlsx` (Table S9, Table S23) |
| DF02 | Mori M. *et al.* (2021) From coarse to fine: the absolute *Escherichia coli* proteome under diverse growth conditions. *Mol. Syst. Biol.* 17:e9536. [doi:10.15252/msb.20209536](https://doi.org/10.15252/msb.20209536) | absolute protein mass fraction | `2.3MOESM8_ESM.xlsx` (EV8), `2.6MOESM2_ESM.xlsx` (EV2) |
| DF03 | Wu C. *et al.* (2023) Enzyme expression kinetics by *Escherichia coli* during transition from rich to minimal media depends on proteome reserves. *Nat. Microbiol.* 8:347–359. [doi:10.1038/s41564-022-01310-w](https://doi.org/10.1038/s41564-022-01310-w) | absolute protein mass fraction | `3.1CHEN41564_2022_1310_MOESM3_ESM.xlsx` (Group #1–#3) |
| DF04 | Zhu M., Mori M., Hwa T., Dai X. (2025) Distantly related bacteria share a rigid proteome allocation strategy with flexible enzyme kinetics. *PNAS*. [doi:10.1073/pnas.2427091122](https://doi.org/10.1073/pnas.2427091122) | DIA, absolute protein mass fraction | `4.1. Zhu_SD01.xlsx`, `4.2Zhu_SD02.xlsx` (E. coli sheets) |

The supplementary files are included in `data/raw/` exactly as downloaded, so the analysis can be
rerun without any download. Their SHA-256 checksums are recorded in `results/run_log.json`.
The other files in `data/raw/` (`2.2…`, `2.4…`, `2.5…`) are additional supplementary tables that
are not read by the pipeline.

## Analysis workflow

```
raw supplementary tables ──► 01 ingest & harmonise ──► 02 metadata, scope & QC ──► 03 COG allocation vs growth (primary)
                                                                               ├─► 04 ASCA, protein level (secondary)
                                                                               ├─► 05 prediction, leave-one-dataset-out
                                                                               └─► 06 sensitivity to missing values
                                                         all saved results ──► 07 final figures
```

| step | notebook | what it does |
|---|---|---|
| 1 | `01_ingest_harmonise` | one loader per paper → one samples × b-number matrix of protein **mass fractions**; one COG annotation for all datasets; missingness QC |
| 2 | `02_metadata_scope_qc` | harmonised metadata table; every sample classified by perturbation type; scope = nutrient-quality conditions; design and unit checks |
| 3 | `03_cog_allocation_vs_growth` | **primary result**: COG-category mass fraction vs growth rate (within-dataset slope, permutation p, bootstrap CI, FDR, leave-one-dataset-out stability) |
| 4 | `04_asca_protein_level` | ASCA on log protein abundances: variance partition dataset / growth / interaction / residual, permutation tests, COG-level summary |
| 5 | `05_prediction_validation` | PLS regression of growth rate on COG fractions, nested leave-one-dataset-out cross-validation, mean baseline, permutation null, jackknife coefficients |
| 6 | `06_sensitivity_censoring` | repeats the notebook-03 model under four treatments of left-censored (undetected) values |
| 7 | `07_figures` | redraws all final figures from the saved results in one consistent style (no statistics recomputed) |

Every notebook starts with a markdown cell explaining **why** the step is needed, and each code
cell is preceded by the reasoning behind it. The functions the notebooks call are in
`src/proteome_alloc/` and each has a docstring explaining what it does and why.

## Repository layout

```
README.md               this file
config.yml              every parameter: files, sheets, columns, thresholds, seed, number of workers
environment.yml         software environment (Python and package versions)
run_all.py              single entry point: annotation + notebooks 01–07, in order
data/raw/               supplementary files as downloaded (read-only, see data/raw/README.md)
data/annotation/        cached UniProt + COG tables and provenance.json (URL, date, release, SHA-256)
data/processed/         harmonised matrix and metadata (generated by notebooks 01–02)
metadata/               curated sample tables, COG category names, documented scope overrides
src/proteome_alloc/     all functions (config, ingest, annotation, scope, allocation, asca,
                        combined, predict, proteomap, plots)
notebooks/              seven notebooks that call src/ and explain each step
results/figures/        figures written by notebooks 01–06; final/ = figures from notebook 07
results/tables/         result tables (TSV), prefixed with the notebook number
results/run_log.json    per step: parameters, input checksums, package versions, outputs
```

Raw data, code and generated results are kept in separate folders.

## How to run

**Requirements:** Python 3.13 with the packages pinned in `environment.yml`
(numpy, pandas, scipy, scikit-learn, matplotlib, openpyxl, PyYAML, joblib, requests,
plus JupyterLab/nbconvert to run the notebooks). Tested on Windows 11, on VS code software

```bash
git clone https://github.com/AlfArlin/E.coli.git
cd E.coli
conda env create -f environment.yml
conda activate proteome-allocation

python run_all.py              # runs everything, about 2 minutes on a laptop
python run_all.py --from 3     # restart from notebook 03
```

`run_all.py` executes each notebook top to bottom in a fresh kernel and saves it in place, so the
committed notebooks always show the outputs of a clean, complete run. Alternatively, open the
notebooks in JupyterLab and run them one by one **in numerical order** (each one reads the files
written by the previous ones).

No manual step and no internet connection are needed: the data and the annotation cache are in
the repository. `python run_all.py --refresh-annotation` re-downloads the UniProt and COG files
(this changes the annotation if the databases have been updated since the cached release).

## Conclusions and what I have learned

**Merging processed datasets is the main difficulty.** The four studies use different units
(Schmidt reports copies per cell, the others mass fractions), different instruments and
different detection depths. Only 1,188 of 4,449 proteins are quantified in every
in-scope sample. Growth rates and sample metadata had to be collected by hand
from each paper, because they are not in a common format.

**Why the prediction model is not yet reliable.** On average, the PLS model predicts growth
rate in a held-out study better than the mean (RMSE 0.36 vs 0.53 h⁻¹, permutation
p = 0.005), but it does not transfer to every study. The significant dataset × growth-rate interaction in
ASCA (p = 0.013) shows that the studies do not fully agree on how the proteome responds,
which limits transfer between them.

**Preliminary conclusion.** The categories associated with growth rate follow the known
bacterial growth laws. Translation and ribosomal proteins (COG J) increase most strongly, together
with chaperones/turnover (O), cell envelope (M) and transcription (K). Amino acid metabolism
(E), energy production (C), signal transduction (T), inorganic ion and carbohydrate
transport (P, G) and motility (N) decrease. All these slopes keep their sign when any one
dataset is left out, and the 14 significant categories are the same under all four
treatments of missing values (notebook 06). At protein level, amino acid biosynthesis
enzymes (*leuB*, *ilvB*, *ilvC*, *argD*, *argG*, *serA*, *gdhA*) and uptake systems
(*dppA*, *gltI*, *argT*) fall most with growth rate, while *ompC*, *rapA*, *deaD* and *fis*
rise. These are associations across 28 conditions, not causal effects, and they need
confirmation in further independent studies.
## Reproducibility

* **All parameters in `config.yml`**; none are hidden in notebooks or functions.
* **Initially this work was written in 1 notebook** eventually they were split to keep a logical order on github
* **On seed** (`seed` in `config.yml`). Each stochastic step drives from that seed.
* **Parallel computing.** Bootstrap resamples and permutation tests run in parallel with
  `joblib` (`n_jobs` in `config.yml`).
* **Run record.** 


## Design decisions 

* **Same unit in all datasets.** Schmidt reports copies per cell. It is converted with
  `N·MW / Σ(N·MW)` (fraction of measured protein mass), the quantity the other studies
  report. No cell-volume or dry-weight constant is needed.
* **One COG annotation.** b-number → COG id (UniProt eggNOG cross-reference) → category
  letter (COG 2020). A protein with several categories gives an equal share of its mass
  to each; proteins without a COG are kept as group `-`.
* **No imputation in the primary analysis.** Category mass fractions are sums, and an
  undetected protein contributes zero mass. 
* **Scope: nutrient quality.** Batch culture, 37 °C, no stress, antibiotic, stationary phase
  or special procedure. Rules are regular expressions in `config.yml, and saved in
  `results/tables/02_scope_decisions.tsv`.
  (no imputation). Preprocessing order: filter → log10 → average replicates → centre.
* **Validation without leakage.** Leave-one-dataset-out, nested (the number of PLS components
  is chosen on training data only); scaling is inside the scikit-learn pipeline.
## Known limitations

* Schmidt 2016 contributes two NCM3722 conditions in scope, i.e. a single within-dataset contrast.
* Growth rates for Wu 2023 and Zhu 2025 were entered by hand in `metadata/wu_2023_samples.tsv`
  and `metadata/zhu_2025_samples.tsv`; the `Growth_Rate_Source` column gives their origin.
  Mori 2021 growth rates are computed from the published doubling times, Schmidt 2016 from Table S23.
* The scope classification of Mori 2021 is keyword-based (reviewed decisions are in
  `metadata/scope_overrides.tsv`).
* Category fractions sum to 1, so a falling category can be passive dilution by a growing one
  (notebook 03, section 4).
* Internal validation on four studies of one strain background; no claim is made beyond
  nutrient-limited growth of *E. coli* K-12.


  ## Use of AI tools
Transparency statement: this project was developed with the help of an AI assistant,
Gemini from google tools during the course (September–October 2026).

**What the AI assistant was used for**

* Discussing analysis options and their statistical justification (e.g. how to treat missing
  values, ASCA, leave-one-dataset-out validation, avoiding data leakage).
* Drafting and restructuring code: the dataset loaders, the functions in `src/proteome_alloc/`,
  and the split of the original single notebook into the numbered notebooks were written with
  the assistant and revised over several iterations.
* Debugging errors.
* Documentation: docstrings, comments, the markdown explanations in the notebooks and this
  README, including a final review of the repository against the course checklist
  (18 October 2026), in which the pipeline was rerun on a second operating system.

**What I decided and did myself**

* The research question, the choice of the four datasets and the scope (nutrient-quality
  conditions only).
* Discussing analysis options during course meetings and their statistical justification ( how to treat missing
  values, ASCA, leave-one-dataset-out validation, avoiding data leakage).
* The key analysis decisions: one COG annotation for all datasets, no batch-correction
  algorithm, growth condition as the statistical unit.
* Collecting the supplementary files and curating the sample metadata in `metadata/`.
* Running the pipeline, inspecting every output and figure, and interpreting the results.

**How the AI-assisted output was checked**

* Every notebook runs top to bottom from the raw files, and the code stops on inconsistent
  sample IDs, missing metadata or non-orthogonal effect matrices.
* Converted values were checked against the original publications (e.g. the Schmidt
  copies-per-cell conversion and the EF-Tu / HU cross-check in notebook 01).
* The full pipeline was reproduced on a second operating system with identical p-values.

I reviewed all code and text suggested by the assistant and take responsibility for the final
content, including any remaining errors.
