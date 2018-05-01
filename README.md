# Chronic myeloid leukaemia gene expression app

## Introduction 
Chronic-phase Chronic myeloid leukaemia (CP-CML) is a blood cancer that affects the blood and bone marrow. In CML, it produces high white blood cells counts that interfere with normal blood cell production. CML is characterised by translocation of chromosome 9 and 22 which resulting as a fusion gene BCR-ABL1. Although imatinib has revolutionised the treatment of CML, only 60% of patients will have an excellent response. Of the remainder, some will progress to more aggressive, fatal blast crisis phase (BC).


## Aim
The aim and motivation of this project was to identify gene expression that were associated with BC-CML patients compared to the CP-CML patient samples, based on gene expression data generated by the next-generation sequencing. Identification of differentially gene expression could aid understanding the biology of BC-CML, potentially leading to new drug target and novel treatment discovery.

## Data source
Public repository database, NCBI sequence read archive ([SRA](https://www.ncbi.nlm.nih.gov/sra)), was used to identify a BC-CML dataset with an accession number [SRP028528](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP028528). There were 9 BC-CML and 8 CP-CML patient samples in this dataset. 

## R packages installation requirement
* shinydashboard
* DT
* shiny
* ggplot2
* gplots
* tidyverse
* heatmaply
* RColorBrewer
* plotly
* networkD3
* igraph
* limma (Bioconductor)

## Example of installing above R packages
```r
# install required packages
install.packages(c("shinydashboard", "DT", "shiny", "ggplot2", "gplots", 
                   "tidyverse", "heatmaply", "RColorBrewer", "plotly", 
                   "networkD3", "igraph"))

# install required bioconductor package (https://bioconductor.org)
source("http://bioconductor.org/biocLite.R") 
biocLite("limma")
```

## Instruction of running the visualisation product (R shiny dashboard)
1. Download the *R* source code here [(https://github.com/chungkok/CML_GEP_App)](https://github.com/chungkok/CML_GEP_App)
2. The user requires to unzip the folder and set to the targeted directory in Rstudio or R. 
3. Prior to run the app, the user must have all the R packages mentioned above installed. 
4. Then the user could type `shiny::runApp()` to start the app.

## Things to note
Best view with wide screen and with at least 8GB RAM pc or mac (23776 genes to be explored). Most figures or plots can be zoom in and out, download as png format (except for 3D scatterplot), pan, and rotation (for multidimensional plot). The toolbar will be appeared when hover to the image/plot on the top right corner (small icons).


