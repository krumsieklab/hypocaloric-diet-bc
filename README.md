# Hypocaloric Diets and Breast Cancer risk

Schweickart *et al.*, "Metabolomic signatures of hypocaloric dietary interventions associate with breast cancer risk in the Nurses’ Health Study II", unpublished.


## File overview


| File                                       | Description                                                                               |
|--------------------------------------------|-------------------------------------------------------------------------------------------|
| `internal_functions.R`                     | Helper functions (e.g. regression utilities) used in other scripts                       |
| `01_Preprocess_Controlled_Trials_data.R`   | Harmonizes and preprocesses data from MSK & OSU controlled trials                        |
| `02_get_differential_metabolites.R`        | Identifies diet-dependent differential metabolites                                       |
| `03_Preprocess_NHS_data.R`                 | Preprocesses the NHSII data                                                              |
| `04_load_save_data.R`                      | Calculates diet scores and merges with NHSII metadata                                    |
| `05_run_mwas.R`                            | Performs metabolite–breast cancer risk MWAS                                              |
| `06_calculate_metabolite_concordance.R`    | Computes concordance of metabolite effects                                               |
| `Schweickart_2025_Tables.Rmd`              | Generates summary tables from processed data                                             |
| `Schweickart_2025_Figures.Rmd`             | Produces Figures 2-4.                                      |

To knit the R Markdown files:

- `rmarkdown::render("Schweickart_2025_Tables.Rmd")`  
- `rmarkdown::render("Schweickart_2025_Figures.Rmd")`

## Note on dataset 

This repository does not contain the input data files for these scripts. Contact https://nurseshealthstudy.org/contact for NHSII data, and the corresponding authors for OSU and MSK data.
