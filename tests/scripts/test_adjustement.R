# Set seed for reproducibility
set.seed(42)

# Function for Holm-Bonferroni adjustment (like your C++ code)
adjusted_holm <- function(p_values) {
  m <- length(p_values)
  indexed <- data.frame(p = p_values, index = seq_along(p_values))
  indexed <- indexed[order(indexed$p), ]

  adjusted <- numeric(m)
  prev <- 0
  for (i in seq_len(m)) {
    raw <- min((m - i + 1) * indexed$p[i], 1)
    adjusted[i] <- max(prev, raw)  # enforce monotonicity
    prev <- adjusted[i]
  }

  # Reorder to original
  reordered <- numeric(m)
  reordered[indexed$index] <- adjusted
  return(reordered)
}

# Define parameters
intervals <- list(
  c(0.1, 1),
  c(0.01, 0.1),
  c(0.001, 0.01)
)
ks <- c(2, 3, 4)

# Store results
results <- list()
counter <- 1

# Simulate and adjust
for (intv in intervals) {
  for (k in ks) {
    n <- 100000 * k
    p_vals <- runif(n, min = intv[1], max = intv[2])
    p_adj <- adjusted_holm(p_vals)
    
    results[[counter]] <- list(
      p = p_vals,
      p_adj = p_adj,
      label = paste0("Interval [", intv[1], ", ", intv[2], "] - k=", k)
    )
    
    counter <- counter + 1
  }
}

# Plotting
library(ggplot2)

plot_list <- list()

for (i in seq_along(results)) {
  df <- data.frame(
    p_value = c(results[[i]]$p, results[[i]]$p_adj),
    type = rep(c("Original", "Adjusted"), each = length(results[[i]]$p)),
    label = results[[i]]$label
  )
  
  p <- ggplot(df, aes(x = p_value, fill = type)) +
    geom_histogram(bins = 100, alpha = 0.6, position = "identity") +
    facet_wrap(~type, ncol = 1) +
    ggtitle(results[[i]]$label) +
    theme_minimal()
  
  plot_list[[i]] <- p
}

# Display all plots (you can use patchwork or gridExtra to arrange them)
library(patchwork)
wrap_plots(plot_list, ncol = 1)
