# --- Load libraries ---
args <- commandArgs(trailingOnly = TRUE)
library(readr)
library(dplyr)

# --- Load input files from command line ---
feature_file <- args[1]
pheno_file <- args[2]

# --- Load features ---
features <- read_tsv(feature_file)
colnames(features)[1] <- "IID"

# --- Load phenotype ---
pheno <- read_tsv(pheno_file, col_types = cols_only(IID = col_character(), PHENO = col_double()))

# --- Merge on sample ID ---
df <- inner_join(features, pheno, by = "IID")

# --- Prepare data for model ---
X <- df %>% select(-IID, -PHENO)
y <- df$PHENO

# --- Add intercept manually (R does this by default in lm) ---
model <- lm(y ~ ., data = X)

# --- Print summary ---
summary(model)

# Run it in terminal like:
# remove the last column header if column empty
# Rscript linear_regression.R ../output/regression/48_51.tsv ../data/quantitative/phenotype.tsv
