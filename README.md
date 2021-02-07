# Analysis of the Relationship Between Essential and Non-Essential Schizosaccharomyces pombe Genes Reveals a Correlation Between mRNA Function and localisation on mRNA Stability

## Project Description:

This project is based on the big data biology assessment completed last stage. The aim is to rework the code used in the previous report to make it neater and reproducible and covert the previous report from a word document to an Rmd file.

Due to improvements in the import code, all of the localisation data (Matsuyama et al. 2006) has been included rather than the subset used before. This resulted in a mild change to the outcome of the statistical tests. These changes warranted slight edits from the original report written in word (see report/"Analysis of the effect of function on mRNA stabilities.docx"). These changes are minor and are just to remove/edit statements and suggestions that are no longer valid/relevant. If there are any statements left that do not make sense just disregard them; the focus of the project is just to convert the word document into an Rmd, not to build upon the previous analysis.

A package called "presetplots" was created in order to compile created functions. These functions aim to improve readability of the code, save custom plot preferences and combine multiple plots into a single function in order to de-clutter the code.

For comparison there is an annotated version of both the old (scripts/"old-master") and new (scripts/"new-master") scripts used to create the Rmd document. The new import and analysis code have also been separated into separate files in order to make them more readable. 

## Technical description
R version 4.0.3 (2020-10-10)

### Packages
-   car_3.0-10
-   carData_3.0-4
-   Rmisc_1.5
-   plyr_1.8.6
-   lattice_0.20-41
-   ggpubr_0.4.0
-   viridis_0.5.1
-   viridisLite_0.3.0
-   patchwork_1.1.1
-   forcats_0.5.0
-   stringr_1.4.0
-   dplyr_1.0.2
-   purrr_0.3.4
-   readr_1.4.0
-   tidyr_1.1.2
-   tibble_3.0.4
-   ggplot2_3.3.3
-   tidyverse_1.3.0
-   hrbrthemes_0.8.0
-   presetplots_0.1.0

The output of `sessionInfo()` can be found in [Session Information](sessioninfo.md)

### Notes
In order to run any of the code you need to download the file "AnGeLiDatabase.txt" and save it to the data-raw folder. This file is too large to send directly.

