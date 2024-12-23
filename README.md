---
title: "Readme"
output: html_document
---

This is the readme file for the article "Geometric microliths as cultural markers: shape/size variation and projectile performance", by A. Cortell-Nicolau, A. Key and A. Palomo.

This folder contains the following subfolders: Data, Figures, Figures_supplementary, Outlines, Results, Scripts, Shapes and Utility_scripts. For a quick check, ideally, only loading and running the script `Reproducible_script.R` within the folder Scripts should work. If you run it with the current number of simulations (300,000) it can take long, but you can try with less simulations (change lines 21-25 in that script) to whatver suits you better for a first check. In detail, each folder contains:

* **Data:** This contains the raw and clean data, as follows:
  * *Exp_data.csv:* This is the data as it was downloaded from the tensile tester, combined with the metric data that was gathered during the experiment. It also includes some utility information (e.g. geometric original ID from the experiment was copied).
  * *Full_data.csv:* This is the same as above, but including two covariates (PC1 and PC2) from the GMM analysis.
  * *Geo_throws_D1.csv:* Metric data gathered for each geometric during day/shot 1.
  * *Geo_throws_D2.csv:* Metric data gathered for each geometric during day/shot 2.
  * *Geo_throws_D3.csv:* Metric data gathered for each geometric during day/shot 3.
  * *Test_1.csv:* Raw data from the tensile tester for each geometric during day/shot 1.
  * *Test_2.csv:* Raw data from the tensile tester for each geometric during day/shot 2. 
  * *Test_3.csv:* Raw data from the tensile tester for each geometric during day/shot 3.
  * *Results_clean.csv:* Clean merged data. In this file, missing shots (e.g. not properly recorded) have been eliminated.
  * *Sel_data.csv:* Final dataset. **This is the dataset used for the analisys**.

* **Figures:** Figures used in the article (only figures produced in R).

* **Figures_supplementary:** Figures with model diagnostics, etc. These contain the different models tried for each response.

* **Outlines:** .jpg outlines for GMM analysis.

* **Results:** Data obtained for the results of the models and used for final reporting and plots.

* **Scripts:** Contains the script `Reproducible_script.R`. As said before, this should be enough to obtain the results reported in the paper.

* **Shapes:** Original spatial shapes to transform to jpg outlines as per the protocol described in Cortell-Nicolau, 2019.

* **Utility_scripts:** Contains additional scripts used for exploratory analyses, data preparation, model comparison and variable assessment.

*References*

Cortell-Nicolau, A. (2019). *Geomeasure: GIS and scripting for measuring morphometric variability*, Lithic Technology, 44(3), 153-165.
