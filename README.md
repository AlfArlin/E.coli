##Arlin Alfaro Arias -MSB1015 - October 
#Proteome allocation protein groups involve ingrowth rate in *E. coli* 

** Research Question.** Which functional protein of the *E. coli* _K-12 _ proteome are associated with growth rate changes across nutrient conditions,
and does that allocation pattern transfer between independent studies?

**Data.** Four published, already-processed absolute proteomics datasets:
**DF01**Schmidt 2016, Mori 2021, Wu 2023 and Zhu 2025. No raw MS data are reprocessed. If you want to check the references go to the last section.

**Approach.**

| step | notebook | what it does |
|---|---|---|
| 1 | `01_ingest_harmonise` | one loader per paper → one samples × b-number matrix of protein **mass fractions**; one COG annotation for all datasets |
| 2 | `02_metadata_scope_qc` | harmonised metadata table; every sample classified by perturbation type; scope = nutrient-quality conditions; design and unit checks |
| 3 | `03_cog_allocation_vs_growth` | **primary result**: COG-category mass fraction vs growth rate (within-dataset slope, permutation p, bootstrap CI, FDR, leave-one-dataset-out stability) |
| 4 | `04_asca_protein_level` | ASCA on log protein abundances: variance partition dataset / growth / residual, permutation test, COG-level summary |
| 5 | `05_prediction_validation` | PLS regression of growth rate on COG fractions, nested leave-one-dataset-out cross-validation, baseline, permutation null, jackknife coefficients |

## Repository layout

```
config.yml              every parameter: files, sheets, columns, thresholds, seed
environment.yml         software environment
run_all.py              single entry point (annotation + five notebooks, in order)
data/raw/               supplementary files as downloaded (not modified, see data/raw/README.md)
data/annotation/        cached UniProt + COG tables and provenance.json
data/processed/         harmonised matrix and metadata (generated)
metadata/               curated sample tables, COG category names, scope overrides
src/proteome_alloc/     all functions (ingest, annotation, scope, allocation, asca, predict, plots, config)
notebooks/              five thin notebooks that call src/ and explain each step
results/figures/        generated figures (PNG)
results/tables/         generated tables (TSV)
results/run_log.json    parameters, input checksums, package versions per step
tests/                  synthetic-data generator and end-to-end smoke test
```

Raw data, code and generated results are kept in separate folders.

## How to run

```bash
conda env create -f environment.yml
conda activate proteome-allocation

#To run all the pipeline complete
python run_all.py               # builds the COG table (first run needs internet) and runs notebooks 01-05
```

`run_all.py` executes the notebooks in place with `jupyter nbconvert`, so the saved
notebooks contain the outputs of a clean top-to-bottom run. The notebooks can also be
opened and run one by one in Jupyter, in numerical order.

The first run downloads two annotation files (UniProt, E. coli K-12, Swiss-Prot; NCBI
COG 2020 definitions) (Riley, 2006; Rudd, 2000).. 
They are cached in `data/annotation/` together with URL, date,
UniProt release and SHA-256. 

This is the connection of uniprot with the notebooks, commit them just once, later runs are then offline and unaffected by
database updates. If the download is blocked, place the files manually (same file names
as in `config.yml`).

## Design decisions:

* **Same unit in all datasets.** DF01 reports copies per cell. It is converted with
  `N·MW / Σ(N·MW)` (fraction of measured protein mass), the quantity the other studies
  report. No cell-volume or dry-weight constant is needed.
* **One COG annotation.** b-number → COG id (UniProt eggNOG cross-reference) → category
  letter (COG 2020). A protein with several categories gives an equal share of its mass
  to each (`annotation.multi_category`); proteins without a COG are kept as group `-`.
* **No imputation in the primary analysis.** Category mass fractions are sums, and an
  undetected protein contributes zero mass, therefore we are not imputing proteins that might
  not be in all datasets due to LLOD.
* **Scope: nutrient quality.** Batch culture, 37 °C, no stress, antibiotic, stationary
  or special procedure. The construction of the metadata is structured in `config.yml.
  Result of this are in results/tables:01_ingest_report.
* **Condition = statistical unit.** Replicate runs of one medium are averaged first.
* **Dataset effect.** Each dataset keeps its own intercept; growth effects are estimated
  from within-dataset contrasts only. No batch-correction algorithm is applied to the data.
* **ASCA protein set.** Proteins quantified in every in-scope sample of every dataset
  (no imputation). Preprocessing order: filter → log10 → average replicates → centre.
* **Validation.** Leave-one-dataset-out, nested; scaling inside the pipeline; mean
  baseline; permutation null; every tuning trial saved.

## Reproducibility

* All parameters in `config.yml`; none in notebooks or functions.
* One master seed; each stochastic step uses a named stream derived from it, so results
  do not depend on the number of parallel workers (`n_jobs`).
* Paths are relative to the repository root.
* `results/run_log.json` records, per step: config checksum, parameters, SHA-256 of every
  input file, Python/OS/package versions, outputs and all drop counts.
* Assertions stop the run when sample IDs do not match, a metadata row is missing, or
  effect matrices are not orthogonal; warnings report samples without a data column and
  profiles that do not sum to ≈ 1.
* `tests/smoke_test.py` runs the whole pipeline twice on synthetic data and checks that
  all result tables are byte-identical and that a planted signal is recovered.

## Data sources

| dataset | files expected in `data/raw/` | original unit |
|---|---|---|
| Schmidt et al. 2016 | `1.1The quantitative and condition-dependent Escherichia coli proteome.xlsx` (Table S9, Table S23) | copies per cell |
| Mori et al. 2021 | `2.3MOESM8_ESM.xlsx` (EV8), `2.6MOESM2_ESM.xlsx` (EV2) | mass fraction |
| Wu et al. 2023 | `3.1CHEN41564_2022_1310_MOESM3_ESM.xlsx` (Group #1–#3) | mass fraction |
| Zhu et al. 2025 | `4.1. Zhu_SD01.xlsx`, `4.2Zhu_SD02.xlsx` | mass fraction |

Add the full citations and download links of the four papers here before submission.

## Known limitations

* Schmidt contributes two NCM3722 conditions, i.e. one within-dataset contrast.
* Growth rates and media for Wu 2023 and Zhu 2025 are in `metadata/*.tsv`; the
  `Growth_Rate_Source` column must cite the table of the paper they were read from.
* The Mori scope classification is keyword-based and must be reviewed once against the paper.
* Category fractions sum to 1: a falling category can be passive dilution (see notebook 03, section 4).
* Internal validation on four studies of one strain; no claim beyond nutrient-limited growth.


##References 

Galperin, M. Y., Wolf, Y. I., Makarova, K. S., et al. (2020). COG database update: focus on microbial diversity, model organisms, and widespread pathogens. Nucleic Acids Research, 49, D274–D281. https://doi.org/10.1093/nar/gkaa1018
