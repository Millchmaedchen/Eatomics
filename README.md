# Eatomics
Eatomics is a web application developed using Shiny. Eatomics enables fast exploration of differential expression and pathway analysis to researchers with limited bioinformatics knowledge.The application will aid in quality control of the quantitative proteomics data, visualization, differential expression and pathway analysis. Highlights of the application are an extensive experimental setup module, the data and report generation feature and the multiple ways to interact and customize the analysis. 

# Getting started
## 1. Install Rstudio
If you do not have installed Rstudio yet, please find the installation for your distribution at https://rstudio.com/products/rstudio/download/ . 
## 2. Run Shiny App in Rstudio
Within R studio please run

```
install.packages("shiny")
install.packages("devtools")
install.packages("pacman")

library(shiny)
library(devtools)
library(pacman)

runUrl("https://github.com/Millchmaedchen/Eatomics/archive/master.zip", subdir = "R")
```
## 3. Input Files
Two demo data files are downloaded when the above runUrl() command is used. You may either find those files in the local folder in Eatomics/Data or you can download them separately from this repository. 
1. Demo_proteinGroups.txt: The proteinGroups.txt (i.e. a tab-separated files) as generated by the quantitative analysis software of raw mass spectrometry data - MaxQuant. The file should contain at least the columns Protein IDs, Majority protein IDs, Gene names, LFQ/iBAQ measurement columns, Reverse, Potential contaminant, Only identified by site. The latter three may be empty.
2. Demo_clinicaldata.txt: The sample description file - a tab separated text file as can be produced with any Office program by saving the spread sheet as .txt. The file shoul contain an m x n matrix file containing the clinical/observational/experimantal details of m number of patients with n number of parameters. 

## 4. Troubleshooting 

- R version < 3.6 will need a lower version of caTools in order to create proper reports from rmarkdown - this might help to fix any caTools related errors.
```
install.packages("https://cran.r-project.org/src/contrib/Archive/caTools/caTools_1.14.tar.gz", repos=NULL, type="source")
```

# References

1. Krug, Karsten, et al. "A curated resource for phosphosite-specific signature analysis." Molecular & cellular proteomics 18.3 (2019): 576-593. - describes ssGSEA 2.0, which is re-used in this repository
2. Ritchie, Matthew E., et al. "limma powers differential expression analyses for RNA-sequencing and microarray studies." Nucleic acids research 43.7 (2015): e47-e47. - describes limma, which is needed for the core functionality of differential expression analysis
