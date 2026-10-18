# Raw data (not modified by the pipeline)

The supplementary files of the four publications, exactly as downloaded from the journals.
The pipeline only reads them; file names and sheet names are set in `config.yml`.
The SHA-256 of every file used is written to `results/run_log.json` when notebook 01 runs,
so anyone can check they are using identical files.

| dataset | publication (DOI) | file | sheets used |
|---|---|---|---|
| DF01 Schmidt 2016 | [10.1038/nbt.3418](https://doi.org/10.1038/nbt.3418) | `1.1The quantitative and condition-dependent Escherichia coli proteome.xlsx` | `Table S9` (copies per cell), `Table S23` (growth rates) |
| DF02 Mori 2021 | [10.15252/msb.20209536](https://doi.org/10.15252/msb.20209536) | `2.3MOESM8_ESM.xlsx` | `EV8-AbsoluteMassFractions-1` |
| DF02 Mori 2021 | | `2.6MOESM2_ESM.xlsx` | `EV2-Samples-1` (sample annotation, doubling times) |
| DF03 Wu 2023 | [10.1038/s41564-022-01310-w](https://doi.org/10.1038/s41564-022-01310-w) | `3.1CHEN41564_2022_1310_MOESM3_ESM.xlsx` | `Group #1/#2/#3 protein mass fractions` |
| DF04 Zhu 2025 | [10.1073/pnas.2427091122](https://doi.org/10.1073/pnas.2427091122) | `4.1. Zhu_SD01.xlsx`, `4.2Zhu_SD02.xlsx` | E. coli sheets; all sheets are searched for the sample columns in `metadata/zhu_2025_samples.tsv` |

`2.2MOESM6_ESM.xlsx`, `2.4MOESM9_ESM.xlsx` and `2.5MOESM1_ESM.xlsx` are additional supplementary
tables that are not read by the pipeline.
