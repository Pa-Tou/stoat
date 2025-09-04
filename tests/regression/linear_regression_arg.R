# --- Load libraries ---
args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# --- Load input files ---
feature_file <- args[1]
pheno_file <- args[2]

# Read features (force numeric, fix decimal issues)
features <- read_tsv(feature_file, locale = locale(decimal_mark = ".")) %>%
  mutate(across(-1, as.numeric))  # Convert all but first column to numeric

colnames(features)[1] <- "IID"

# Read phenotype
pheno <- read_tsv(pheno_file, col_types = cols_only(IID = col_character(), PHENO = col_double()))

# Merge
df <- inner_join(features, pheno, by = "IID")

# Prepare data
X <- df %>%
  select(-IID, -PHENO) %>%
  select(where(~ !all(is.na(.)) & length(unique(.)) > 1))  # remove all-NA or constant

y <- df$PHENO

df_model <- X
df_model$PHENO <- y

model <- lm(PHENO ~ ., data = df_model)
summary(model)

# Run it in terminal like:
# remove the last column header if column empty
# Rscript linear_regression_arg.R ../../output/regression/1_4.tsv ../../data/quantitative/phenotype.tsv
# Rscript tests/regression/linear_regression_arg.R output_droso/regression/7165_7162.tsv ../lab/droso/data/pangenome_pheno.tsv
