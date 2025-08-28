library(pheatmap)
library(corrplot)     # for correlation heatmap
library(RColorBrewer) # for color palettes

# Read genotype data (X matrix)
X <- read.table("table_test.tsv", header = FALSE)

# Read phenotype data (Y vector)
Y <- read.table("phenotype.tsv", header = FALSE)[,1]

# Check if lengths match
if (nrow(X) != length(Y)) {
  stop("Mismatch: Number of rows in X and length of Y must be equal.")
}

# Combine into one data frame
df <- cbind(X, Y = Y)

# Convert to data frame
df <- as.data.frame(df)

# Run linear regression: phenotype Y ~ all columns in X
model <- lm(Y ~ ., data = df)

# View regression summary
summary(model)

# Convert to matrix (if not already)
X_matrix <- as.matrix(X)

# Plot heatmap
pheatmap(X_matrix,
         cluster_rows = TRUE,      # cluster samples
         cluster_cols = TRUE,      # cluster features
         show_rownames = FALSE,    # hide row names (samples)
         show_colnames = FALSE,    # hide column names (features)
         main = "Heatmap of Genotype Matrix X")

# Compute correlation matrix
cor_matrix <- cor(X_matrix, use = "pairwise.complete.obs")  # handles missing values if any

# Plot correlation heatmap
corrplot(cor_matrix,
         method = "color",         # use colored squares
         type = "upper",           # only upper triangle
         col = colorRampPalette(brewer.pal(8, "RdBu"))(200),  # color gradient
         tl.col = "black",         # text color
         tl.cex = 0.8,             # text size
         addCoef.col = "black",   # add correlation values
         number.cex = 0.7,         # correlation number size
         diag = FALSE,             # hide diagonal
         title = "Correlation Matrix of Predictors",
         mar = c(0,0,2,0))         # margin for title

print(round(cor_matrix, 2))
# Save correlation matrix to file
write.table(round(cor_matrix, 2), file = "correlation_matrix.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
