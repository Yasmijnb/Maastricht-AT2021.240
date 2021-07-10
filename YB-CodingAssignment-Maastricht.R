################################################################################

# Coding assignment for PhD position AT2021.240 in Maastricht

################################################################################

# Install packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RCy3")
BiocManager::install("biomaRt")

################################################################################

# Import packages
library(biomaRt)
library(RCy3)

################################################################################

# Import data

variants <- read.csv("./variant_list.txt", header = FALSE)

################################################################################

# Step 1: Retrieve the associated gene for each variant using BioMart



# Step 2: Create a network in Cytoscape linking the SCPs to the genes