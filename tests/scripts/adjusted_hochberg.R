# Define a vector of raw p-values
pvals <- c(0.0000491, 0.000025)

# Apply Hochberg correction using base R
hochberg_adj <- p.adjust(pvals, method = "hochberg")

# Show results
data.frame(
  raw_p = pvals,
  hochberg_adjusted_p = hochberg_adj
)