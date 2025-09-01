# Define a vector of raw p-values
pvals <- c(0.00000001, 0.1, 0.32, 0.00002, 0.234, 0.5)

# Apply Hochberg correction using base R
hochberg_adj <- p.adjust(pvals, method = "hochberg")

# Show results
data.frame(
  raw_p = pvals,
  hochberg_adjusted_p = hochberg_adj
)