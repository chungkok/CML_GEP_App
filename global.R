#global.R

# the following packages are required before commerce the program
# install.packages(c("shinydashboard", "DT","shiny", "ggplot2", "gplots", "tidyverse", "heatmaply",
#               "RColorBrewer", "plotly", "networkD3", "igraph"))
#   source("http://bioconductor.org/biocLite.R")
#   biocLite("limma")

## Load required packages
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(gplots)
library(limma)
library(heatmaply)
library(RColorBrewer)
library(tidyverse)
library(plotly)
library(reshape2)
library(networkD3)
library(igraph)

# load sample information
targets <- read.csv(file="./data/targets.csv")
targets$sample_name <- gsub("-", ".", targets$sample_name)

# load transformed counts
logcounts <- read.csv("./data/log transformed counts cpm.csv", stringsAsFactors = FALSE)
logcounts_melt <- melt(logcounts, id.vars="Symbol", variable.name="sample_name")

# logcounts for heatmap and MDS
logcounts_heatmap <- logcounts
rownames(logcounts_heatmap) <- logcounts[,1]
logcounts_heatmap <- logcounts_heatmap[,-1]
var_genes <- apply(logcounts_heatmap, 1, var)

# read in dataset
data <- read.csv(file="./data/topTable BC1 vs BC2 limma voom.csv", stringsAsFactors = FALSE)
data1 <- read.csv(file="./data/topTable BC1 vs CP limma voom.csv", stringsAsFactors = FALSE)
data2 <- read.csv(file="./data/topTable BC2 vs CP limma voom.csv", stringsAsFactors = FALSE)
data3 <- read.csv(file="./data/topTable BC vs CP limma voom.csv", stringsAsFactors = FALSE)
hi2 <- read.delim(file="./PPI_database/combined_PPI.tsv", stringsAsFactors = FALSE)


# data for stripchart
ccat <- merge(logcounts_melt, targets, by.x="sample_name", by.y="sample_name")
ccat <- ccat %>%
  mutate(tooltip = paste(sample_name, "Group: ", group2, "Expression :", signif(value,digits=2), sep = "\n"))

