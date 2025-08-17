set.seed(123)

simulate_pvalues <- function(n = 100, reps = 1000, beta1 = 0, sigma = 1) {
  pvals <- numeric(reps)
  
  for (i in 1:reps) {
    x <- rnorm(n, mean = 0, sd = 1)
    y <- beta1 * x + rnorm(n, mean = 0, sd = sigma)
    
    model <- lm(y ~ x)
    pvals[i] <- summary(model)$coefficients[2, 4]  # p-value for slope
  }
  
  return(pvals)
}

# Simulate under the null (beta1 = 0)
p_null <- simulate_pvalues(beta1 = 0)

# Simulate under the alternative (beta1 = 0.5)
p_alt <- simulate_pvalues(beta1 = 0.5)

# Function to print percentages
print_sig_rates <- function(pvals) {
  cat("p < 0.01:", mean(pvals < 0.01) * 100, "%\n")
  cat("p < 0.05:", mean(pvals < 0.05) * 100, "%\n")
}

cat("Null (beta1 = 0):\n")
print_sig_rates(p_null)

cat("\nAlternative (beta1 = 0.5):\n")
print_sig_rates(p_alt)
