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

### Software Requirements
  
  This code was created with R version 4.2.0 and Rstudio Version 2023.06.2+561 and tested on macOS (Sonoma 14.2.1) with a 2.3 GHz Quad-Core Intel Core i5 CPU.
  
  
### Package Requirements
  
  The following R packages are required to run the above scripts:


  * [maplet (1.1.2)](https://github.com/krumsieklab/maplet)
  * [glue (1.8.0)](https://cran.r-project.org/web/packages/glue/)
  * [dplyr (1.1.4)](https://cran.r-project.org/web/packages/dplyr/)
  * [ggplot2 (4.0.2)](https://cran.r-project.org/web/packages/ggplot2/)
  * [tidyverse (2.0.0)](https://cran.r-project.org/web/packages/tidyverse/)
  * [openxlsx (4.2.8)](https://cran.r-project.org/web/packages/openxlsx)
  * [purrr (1.0.2)](https://cran.r-project.org/web/packages/purrr/)
  * [cowplot (1.2.0)](https://cran.r-project.org/web/packages/cowplot/)
  * [viridis (0.6.3)](https://cran.r-project.org/web/packages/viridis/)
  * [reshape (0.8.10)](https://cran.r-project.org/web/packages/reshape/)
  * [sva (3.46.0)](https://www.bioconductor.org/packages/release/bioc/html/sva.html)
  * [gt (1.0.0)](https://cran.r-project.org/web/packages/gt/)
  * [gtsummary (2.4.0)](https://cran.r-project.org/web/packages/gtsummary/)
  * [ggrepel (0.9.6)](https://cloud.r-project.org/web/packages/ggrepel)
  * [viridisLite (0.4.2)](https://cloud.r-project.org/web/packages/viridisLite/)
  * [data.table (1.14.8)](https://cloud.r-project.org/web/packages/data.table/)
  * [survival (3.8-3)](https://cloud.r-project.org/web/packages/survival/)
  * [lubridate (1.9.4)](https://cloud.r-project.org/web/packages/lubridate/)
  * [forcats (1.0.0)](https://cloud.r-project.org/web/packages/forcats/)
  * [readr (2.1.5)](https://cloud.r-project.org/web/packages/readr/)
  * [tidyr (1.3.1)](https://cloud.r-project.org/web/packages/tidyr/)
  * [tibble (3.2.1)](https://cloud.r-project.org/web/packages/tibble/)
  * [sas7bdat (0.8)](https://cloud.r-project.org/web/packages/sas7bdat/)
  * [chanmetab (0.1.14)](http://dk.archive.ubuntu.com/bioconductor/packages/oldstats/workflows/chanmetab/)
  * [ggpubr (0.6.1)](https://cloud.r-project.org/web/packages/ggpubr/)
  * [caret (7.0-1)](https://cloud.r-project.org/web/packages/caret/)
  * [lattice (0.22-7)](https://cloud.r-project.org/web/packages/lattice/)
  * [glmnet (4.1-10)](https://cloud.r-project.org/web/packages/glmnet/)
  * [Matrix (1.6-4)](https://cloud.r-project.org/web/packages/Matrix/)
  * [stringr (1.5.1)](https://cloud.r-project.org/web/packages/stringr/)
  * [BiocParallel (1.32.6)](https://www.bioconductor.org/packages/release/bioc/html/BiocParallel.html)
  * [genefilter (1.80.3)](https://www.bioconductor.org/packages/release/bioc/html/genefilter.html)
  * [mgcv (1.9-3)](https://cloud.r-project.org/web/packages/mgcv/)
  * [nlme (3.1-168)](https://cloud.r-project.org/web/packages/nlme/)
  * [readxl (1.4.5)](https://cloud.r-project.org/web/packages/readxl/)
  * [SummarizedExperiment (1.28.0)](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
  * [Biobase (2.58.0)](https://www.bioconductor.org/packages/devel/bioc/html/Biobase.html)
  * [GenomicRanges (1.50.2)](https://www.bioconductor.org/packages/devel/bioc/html/GenomicRanges.html)
  * [GenomeInfoDb (1.34.9)](https://www.bioconductor.org/packages/devel/bioc/html/GenomeInfoDb.html)
  * [IRanges (2.32.0)](https://www.bioconductor.org/packages/devel/bioc/html/IRanges.html)
  * [S4Vectors (0.36.2)](https://www.bioconductor.org/packages/devel/bioc/html/S4Vectors.html)
  * [BiocGenerics (0.44.0)](https://www.bioconductor.org/packages/devel/bioc/html/BiocGenerics.html)
  * [MatrixGenerics (1.10.0)](https://www.bioconductor.org/packages/devel/bioc/html/MatrixGenerics.html)
  

<br>


## Note on dataset 

This repository does not contain the input data files for these scripts. Contact https://nurseshealthstudy.org/contact for NHSII data, and the corresponding authors for OSU and MSK data.
