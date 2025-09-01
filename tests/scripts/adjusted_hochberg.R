# Define a vector of raw p-values
pvals <- c(0.02, 0.15, 0.03, 0.001, 0.25, 0.05)

# Apply Hochberg correction using base R
hochberg_adj <- p.adjust(pvals, method = "hochberg")

# Show results
data.frame(
  raw_p = pvals,
  hochberg_adjusted_p = hochberg_adj
)