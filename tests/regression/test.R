# install.packages("RcppEigen")   # if not installed
# install.packages("RcppGSL")     # if not installed

library(Rcpp)
library(RcppEigen)
library(RcppGSL)
library(ggplot2)
library(reshape2)
library(glmnet)

# Load C++ implementations
sourceCpp("/home/mbagarre/Bureau/stoat/tests/regression/linear_regression_maths.cpp")
sourceCpp("/home/mbagarre/Bureau/stoat/tests/regression/linear_regression_stoat.cpp")
sourceCpp("/home/mbagarre/Bureau/stoat/tests/regression/linear_regression_rvtest.cpp")

# Simulation parameters
set.seed(123)
n_sims <- 1
vals <- c(0, 0.5, 1)
p_threshold <- 0.01

# Storage for simulation results
results <- data.frame(
  P_R = numeric(n_sims),
  P_Maths = numeric(n_sims),
  P_Stoat = numeric(n_sims),
  P_RVTest = numeric(n_sims),
  Diff_Maths = numeric(n_sims),
  Diff_Stoat = numeric(n_sims),
  Diff_RVTest = numeric(n_sims),
  N = integer(n_sims),
  K = integer(n_sims)
)

# Create storage for null results
results_null <- data.frame(
  P_R = numeric(n_sims),
  P_Maths = numeric(n_sims),
  P_Stoat = numeric(n_sims),
  P_RVTest = numeric(n_sims),
  Diff_Maths = numeric(n_sims),
  Diff_Stoat = numeric(n_sims),
  Diff_RVTest = numeric(n_sims)
)

# Generate one constrained SNP row (sum of SNPs = 1)
gen_row <- function(num_snps) {
  repeat {
    row <- sample(vals, num_snps, replace = TRUE)
    if (sum(row) == 1) return(row)
  }
}

# ----------------------------------------- SIMULATION ALL -----------------------------------------

for (i in seq_len(n_sims)) {

  # Random sample sizes
  n <- sample(10:100, 1)
  k <- 4 # total predictors incl. intercept [default : 2 path] + intercept
  
  results$N[i] <- n
  results$K[i] <- k
  
  # Step 1a: Generate matrix X (SNP data with sum=1 per row)
  snp_matrix <- matrix(
    t(replicate(n, gen_row(k-1))),
    nrow = n,
    ncol = k-1
  )
  
  # Step 1b: sanity check
  if(any(rowSums(snp_matrix) != 1)) stop("Some rows do not sum to 1!")
  
  print(snp_matrix)
  
  X <- cbind(1, snp_matrix) # Add intercept column
  intercept <- 1
  
  # ---- Step 2a: Significant regression ----
  beta_true <- c(runif(1, -0.5, 0.5), rep(0, ncol(snp_matrix) - 1))
  Y <- intercept + snp_matrix %*% beta_true + rnorm(n, sd = 0.5)
  Y <- as.numeric(Y)
  
  # non-significant regression
  beta_null <- rep(0, ncol(snp_matrix))   # length 2
  Y_null <- intercept + snp_matrix %*% beta_null + rnorm(n, sd = 0.5)
  Y_null <- as.numeric(Y_null)
  
  # ----------------------------------------- SIGNIFICATIVE -----------------------------------------
  
  # Step 3a: Linear regression in R
  df <- data.frame(Y = Y, snp_matrix)
  fit_r <- lm(Y ~ ., data = df)
  r_summary <- summary(fit_r)
  coefs <- coef(r_summary)
  
  # Store the coefficient and p-value of the first predictor
  results$Coef_R[i] <- coef(fit_r)[2]
  results$P_R[i] <- coefs[-1, "Pr(>|t|)"][1] # Exclude intercept p-value
  
  # exclude last column
  X_maths <- snp_matrix[, -ncol(snp_matrix), drop = FALSE]
  
  # Step 4: Linear regression using C++ maths implementation
  res_maths <- cpp_linear_regression_maths(X_maths, Y)
  cpp_coef_maths <- unlist(res_maths$coefficients)
  cpp_pvals_maths <- unlist(res_maths$p_values)
  
  results$P_Maths[i]   <- cpp_pvals_maths[2]
  results$Coef_Maths[i] <- cpp_coef_maths[2]
  
  results$Diff_P_Maths[i]   <- results$P_Maths[i]   - results$P_R[i]
  results$Diff_Coef_Maths[i] <- results$Coef_Maths[i] - results$Coef_R[i]
  
  # Step 5: Linear regression using C++ stoat implementation
  res_stoat <- cpp_linear_regression_stoat(X_maths, Y)
  cpp_coef_stoat <- unlist(res_stoat$coefficients)
  cpp_pvals_stoat <- unlist(res_stoat$p_values)
  
  results$P_Stoat[i]   <- cpp_pvals_stoat[2]
  results$Coef_Stoat[i] <- cpp_coef_stoat[2]
  
  results$Diff_P_Stoat[i]   <- results$P_Stoat[i]   - results$P_R[i]
  results$Diff_Coef_Stoat[i] <- results$Coef_Stoat[i] - results$Coef_R[i]
  
  # Step 6: Linear regression using C++ rvtest implementation
  res_rvtest <- cpp_linear_regression_rvtest(X_maths, Y)
  cpp_coef_rvtest <- unlist(res_rvtest$coefficients)
  cpp_pvals_rvtest <- unlist(res_rvtest$p_values)
  
  results$P_RVTest[i]   <- cpp_pvals_rvtest[2]
  results$Coef_RVTest[i] <- cpp_coef_rvtest[2]
  
  results$Diff_P_RVTest[i]   <- results$P_RVTest[i]   - results$P_R[i]
  results$Diff_Coef_RVTest[i] <- results$Coef_RVTest[i] - results$Coef_R[i]
  
  # ----------------------------------------- NO SIGNIFICATIVE -----------------------------------------
  
  df_null <- data.frame(Y = Y_null, snp_matrix)
  fit_r <- lm(Y ~ ., data = df_null)
  r_summary <- summary(fit_r)
  results_null$P_R[i] <- coef(r_summary)[2, "Pr(>|t|)"]
  results_null$Coef_R[i] <- coef(fit_r)[2]
  
  # C++ Maths
  res_maths <- cpp_linear_regression_maths(snp_matrix[, -ncol(snp_matrix), drop=FALSE], Y_null)
  results_null$P_Maths[i] <- res_maths$p_values[2]
  results_null$Coef_Maths[i] <- res_maths$coefficients[2]
  results_null$Diff_Maths[i] <- results_null$P_Maths[i] - results_null$P_R[i]
  
  # C++ Stoat
  res_stoat <- cpp_linear_regression_stoat(snp_matrix[, -ncol(snp_matrix), drop=FALSE], Y_null)
  results_null$P_Stoat[i] <- res_stoat$p_values[2]
  results_null$Coef_Stoat[i] <- res_stoat$coefficients[2]
  results_null$Diff_Stoat[i] <- results_null$P_Stoat[i] - results_null$P_R[i]
  
  # C++ RVTest
  res_rvtest <- cpp_linear_regression_rvtest(snp_matrix[, -ncol(snp_matrix), drop=FALSE], Y_null)
  results_null$P_RVTest[i] <- res_rvtest$p_values[2]
  results_null$Coef_RVTest[i] <- res_rvtest$coefficients[2]
  results_null$Diff_RVTest[i] <- results_null$P_RVTest[i] - results_null$P_R[i]
  
  results_null$Diff_P_Maths[i]   <- results_null$P_Maths[i]   - results_null$P_R[i]
  results_null$Diff_P_Stoat[i]   <- results_null$P_Stoat[i]   - results_null$P_R[i]
  results_null$Diff_P_RVTest[i]  <- results_null$P_RVTest[i]  - results_null$P_R[i]
  
  results_null$Diff_Coef_Maths[i]  <- results_null$Coef_Maths[i]  - results_null$Coef_R[i]
  results_null$Diff_Coef_Stoat[i]  <- results_null$Coef_Stoat[i]  - results_null$Coef_R[i]
  results_null$Diff_Coef_RVTest[i] <- results_null$Coef_RVTest[i] - results_null$Coef_R[i]
}

# -------------------------------------------------------- PLOT SIGNIFICANCE ----------------------------------------------------

# ---- PLOT p-value distributions ----
df_pvals_long <- melt(
  data.frame(
    P_R      = results$P_R,
    P_Maths  = results$P_Maths,
    P_Stoat  = results$P_Stoat,
    P_RVTest = results$P_RVTest
  ),
  variable.name = "Method",
  value.name = "PValue"
)

ggplot(df_pvals_long, aes(x = PValue, fill = Method)) +
  geom_histogram(alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "P-value Distributions by Method",
       x = "P-value", y = "Frequency")

# ---- Filtered p-value distributions: 1e-1 to 1e-5 ----
df_pvals_long_small <- subset(df_pvals_long, PValue >= 1e-5 & PValue <= 1e-1)

ggplot(df_pvals_long_small, aes(x = PValue, fill = Method)) +
  geom_histogram(alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "P-value Distributions by Method (0 to 1e-5)",
       x = "P-value", y = "Frequency") +
  xlim(1e-5, 1e-1)

# ---- PLOT p-value differences ----
df_pdiff_long <- melt(
  data.frame(
    Diff_P_Maths  = results$Diff_P_Maths,
    Diff_P_Stoat  = results$Diff_P_Stoat,
    Diff_P_RVTest = results$Diff_P_RVTest
  ),
  variable.name = "Method",
  value.name = "PValueDiff"
)

ggplot(df_pdiff_long, aes(x = PValueDiff, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "Difference in P-values vs R lm",
       x = "C++ P-value - R P-value", y = "Frequency")

# ---- PLOT coefficient distributions ----
df_coef_long <- melt(
  data.frame(
    Coef_R      = results$Coef_R,
    Coef_Maths  = results$Coef_Maths,
    Coef_Stoat  = results$Coef_Stoat,
    Coef_RVTest = results$Coef_RVTest
  ),
  variable.name = "Method",
  value.name = "Coefficient"
)

ggplot(df_coef_long, aes(x = Coefficient, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "Coefficient Distributions",
       x = "Coefficient", y = "Frequency")

# ---- PLOT coefficient differences ----
df_diff_long <- melt(
  data.frame(
    Diff_Coef_Maths  = results$Diff_Coef_Maths,
    Diff_Coef_Stoat  = results$Diff_Coef_Stoat,
    Diff_Coef_RVTest = results$Diff_Coef_RVTest
  ),
  variable.name = "Method",
  value.name = "CoefDiff"
)

ggplot(df_diff_long, aes(x = CoefDiff, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "Difference in Coefficients vs R lm",
       x = "C++ Coefficient - R Coefficient", y = "Frequency")

# ---- SIGNIFICANCE proportion ----
p_threshold <- 0.01
sig_props <- sapply(results[, c("P_R", "P_Maths", "P_Stoat", "P_RVTest")],
                    function(p) mean(p < p_threshold, na.rm = TRUE) * 100)

print(sig_props)

# -------------------------------------------------------- PLOT NO SIGNIFICANCE ----------------------------------------------------

# ---- P-value distributions including null ----
df_pvals_long_null <- melt(
  data.frame(
    P_R      = results_null$P_R,
    P_Maths  = results_null$P_Maths,
    P_Stoat  = results_null$P_Stoat,
    P_RVTest = results_null$P_RVTest
  ),
  variable.name = "Method",
  value.name = "PValue"
)

ggplot(df_pvals_long_null, aes(x = PValue, fill = Method)) +
  geom_histogram(alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "P-value Distributions by Method",
       x = "P-value", y = "Frequency")

# ---- PLOT p-value differences ----
df_pdiff_long_null <- melt(
  data.frame(
    Diff_P_Maths  = results_null$Diff_P_Maths,
    Diff_P_Stoat  = results_null$Diff_P_Stoat,
    Diff_P_RVTest = results_null$Diff_P_RVTest
  ),
  variable.name = "Method",
  value.name = "PValueDiff"
)

ggplot(df_pdiff_long_null, aes(x = PValueDiff, fill = Method)) +
  geom_histogram(alpha = 0.5, position = "dodge", bins = 50) +
  theme_bw() +
  labs(title = "Difference in P-values vs R lm",
       x = "C++ P-value - R P-value", y = "Frequency")

# ---- PLOT coefficient distributions ----
df_coef_long_null <- melt(
  data.frame(
    Coef_R      = results_null$Coef_R,
    Coef_Maths  = results_null$Coef_Maths,
    Coef_Stoat  = results_null$Coef_Stoat,
    Coef_RVTest = results_null$Coef_RVTest
  ),
  variable.name = "Method",
  value.name = "Coefficient"
)

ggplot(df_coef_long_null, aes(x = Coefficient, fill = Method)) +
  geom_histogram(alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "Coefficient Distributions",
       x = "Coefficient", y = "Frequency")

# ---- PLOT coefficient differences ----
df_diff_long_null <- melt(
  data.frame(
    Diff_Coef_Maths  = results_null$Diff_Coef_Maths,
    Diff_Coef_Stoat  = results_null$Diff_Coef_Stoat,
    Diff_Coef_RVTest = results_null$Diff_Coef_RVTest
  ),
  variable.name = "Method",
  value.name = "CoefDiff"
)

ggplot(df_diff_long_null, aes(x = CoefDiff, fill = Method)) +
  geom_histogram(alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "Difference in Coefficients vs R lm",
       x = "C++ Coefficient - R Coefficient", y = "Frequency")

# ---- SIGNIFICANCE proportion for non-significant simulations ----
p_threshold <- 0.01

sig_props_null <- sapply(
  results_null[, c("P_R", "P_Maths", "P_Stoat", "P_RVTest")],
  function(p) mean(p < p_threshold, na.rm = TRUE) * 100
)

print(sig_props_null)


# Function to count NA or non-numeric entries in a vector
count_invalid <- function(x) {
  sum(is.na(x) | !is.numeric(x))
}

# For significant results
cat("SIGNIFICANT CASE\n")
cat("OLS: Coef NA/non-numeric =", count_invalid(results$Coef_R), 
    "; P NA/non-numeric =", count_invalid(results$P_R), "\n")
cat("Maths: Coef NA/non-numeric =", count_invalid(results$Coef_Maths), 
    "; P NA/non-numeric =", count_invalid(results$P_Maths), "\n")
cat("Stoat: Coef NA/non-numeric =", count_invalid(results$Coef_Stoat), 
    "; P NA/non-numeric =", count_invalid(results$P_Stoat), "\n")
cat("RVTest: Coef NA/non-numeric =", count_invalid(results$Coef_RVTest), 
    "; P NA/non-numeric =", count_invalid(results$P_RVTest), "\n")

# For non-significant (null) results
cat("\nNULL CASE\n")
cat("OLS: Coef NA/non-numeric =", count_invalid(results_null$Coef_R), 
    "; P NA/non-numeric =", count_invalid(results_null$P_R), "\n")
cat("Maths: Coef NA/non-numeric =", count_invalid(results_null$Coef_Maths), 
    "; P NA/non-numeric =", count_invalid(results_null$P_Maths), "\n")
cat("Stoat: Coef NA/non-numeric =", count_invalid(results_null$Coef_Stoat), 
    "; P NA/non-numeric =", count_invalid(results_null$P_Stoat), "\n")
cat("RVTest: Coef NA/non-numeric =", count_invalid(results_null$Coef_RVTest), 
    "; P NA/non-numeric =", count_invalid(results_null$P_RVTest), "\n")

# ----------------------------------------- SIMULATION COLLINEARITY -----------------------------------------

n_reps   <- 1000

results <- data.frame(
  Rep   = integer(n_reps),
  N     = integer(n_reps),
  K     = integer(n_reps),
  P_R = numeric(n_reps),
  P_Maths = numeric(n_reps),
  P_Stoat = numeric(n_reps),
  P_RVTest = numeric(n_reps),
  Coef_R = numeric(n_reps),
  Coef_Maths = numeric(n_reps),
  Coef_Stoat = numeric(n_reps),
  Coef_RVTest = numeric(n_reps),
  Diff_Maths = numeric(n_reps),
  Diff_Stoat = numeric(n_reps),
  Diff_RVTest = numeric(n_reps),
  Diff_Coef_Maths  <- numeric(n_reps),
  Diff_Coef_Stoat  <- numeric(n_reps),
  Diff_Coef_RVTest <- numeric(n_reps)
)

gen_row <- function(k) {
  # Initialize row with zeros
  row <- rep(0, k)
  # Randomly decide which position(s) will sum to 1
  possible_values <- c(1, 0.5)
  # Randomly pick a combination of values that sum to 1
  repeat {
    row <- sample(c(0, 0.5, 1), size = k, replace = TRUE)
    if (sum(row) == 1) break
  }
  return(row)
}

row_index <- 1   # initialize before loop

for (rep in 1:n_reps) {
  
  # Random sample sizes
  n <- sample(100:1000, 1)
  k <- 5 # total predictors incl. intercept [default : 2 path] + intercept
  number_last_element <- 1  # number of rows to modify at the end
  
  # Step 1: Generate matrix X (SNP data with sum=1 per row)
  X <- matrix(
    t(replicate(n, gen_row(k))),
    nrow = n,
    ncol = k
  )
  
  # Step 2: Add 2 new zero columns
  X <- cbind(X, matrix(0, nrow = n, ncol = 2))
  
  # Step 3: Modify the last 'number_last_element' rows
  if (number_last_element > 0) {
    start_row <- n - number_last_element + 1
    for (i in start_row:n) {
      X[i, ] <- 0                   # reset the whole row to 0
      X[i, (k+1):(k+2)] <- 0.5      # set last 2 columns to 0.5
    }
  }
  
  # Step 3: sanity check
  if(any(rowSums(X) != 1)) stop("Some rows do not sum to 1!")

  # ---- Store meta ----
  results$Rep[row_index] <- rep
  results$N[row_index]   <- n
  results$K[row_index]   <- k
  
  # ---- Generate Y ----
  beta_true <- c(2, rep(0, ncol(X)-1))
  intercept <- 1
  Y <- intercept + X %*% beta_true + rnorm(n, sd = 0.1)

  # ---- R lm ----
  df <- data.frame(Y = Y, X)
  fit_r <- lm(Y ~ ., data = df)
  coefs <- coef(summary(fit_r))
  
  results$P_R[row_index]     <- coefs[-1, "Pr(>|t|)"][1]
  results$Coef_R[row_index]  <- coefs[-1, "Estimate"][1]
  
  dup_cols <- duplicated(as.data.frame(t(X)))
  X_unique <- X[, !dup_cols, drop = FALSE]
  
  # Now proceed with maths version using X_unique instead of X
  X_maths <- X_unique[, -ncol(X_unique), drop = FALSE]
  
  # ---- C++ Maths ----
  res_maths <- cpp_linear_regression_maths(X_maths, Y)
  cpp_pvals_maths <- unlist(res_maths$p_values)
  results$P_Maths[row_index]     <- cpp_pvals_maths[2]
  results$Coef_Maths[row_index]  <- unlist(res_maths$coefficients)[2]
  results$Diff_Maths[row_index]  <- results$P_Maths[row_index] - results$P_R[row_index]
  
  # ---- C++ Stoat ----
  res_stoat <- cpp_linear_regression_stoat(X_maths, Y)
  cpp_pvals_stoat <- unlist(res_stoat$p_values)
  results$P_Stoat[row_index]     <- cpp_pvals_stoat[2]
  results$Coef_Stoat[row_index]  <- unlist(res_stoat$coefficients)[2]
  results$Diff_Stoat[row_index]  <- results$P_Stoat[row_index] - results$P_R[row_index]
  
  # ---- C++ RVTest ----
  res_rvtest <- cpp_linear_regression_rvtest(X_maths, Y)
  cpp_pvals_rvtest <- unlist(res_rvtest$p_values)
  results$P_RVTest[row_index]    <- cpp_pvals_rvtest[2]
  results$Coef_RVTest[row_index] <- unlist(res_rvtest$coefficients)[2]
  results$Diff_RVTest[row_index] <- results$P_RVTest[row_index] - results$P_R[row_index]
  
  # ---- Increment row index ----
  row_index <- row_index + 1
}

# ----------------------------------------- SIMULATION -----------------------------------------

# ---- Filtered p-value distributions: 1e-5 to 0.01 ----
df_pvals_long_small <- subset(df_pvals_long, PValue >= 1e-5 & PValue <= 0.01)

ggplot(df_pvals_long_small, aes(x = PValue, fill = Method)) +
  geom_histogram(alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "P-value Distributions by Method (1e-5 to 0.01)",
       x = "P-value", y = "Frequency") +
  xlim(1e-5, 0.01)

# ---- PLOT p-value distributions ----
df_pvals_long <- melt(
  data.frame(
    P_R      = results$P_R,
    P_Maths  = results$P_Maths,
    P_Stoat  = results$P_Stoat,
    P_RVTest = results$P_RVTest
  ),
  variable.name = "Method",
  value.name = "PValue"
)

ggplot(df_pvals_long, aes(x = PValue, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "P-value Distributions by Method",
       x = "P-value", y = "Frequency")

# ---- PLOT p-value differences ----
df_pdiff_long <- melt(
  data.frame(
    Diff_P_Maths  = results$Diff_Maths,
    Diff_P_Stoat  = results$Diff_Stoat,
    Diff_P_RVTest = results$Diff_RVTest
  ),
  variable.name = "Method",
  value.name = "PValueDiff"
)

ggplot(df_pdiff_long, aes(x = PValueDiff, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "Difference in P-values vs R lm",
       x = "C++ P-value - R P-value", y = "Frequency")

# ---- PLOT coefficient distributions ----
df_coef_long <- melt(
  data.frame(
    Coef_R      = results$Coef_R,
    Coef_Maths  = results$Coef_Maths,
    Coef_Stoat  = results$Coef_Stoat,
    Coef_RVTest = results$Coef_RVTest
  ),
  variable.name = "Method",
  value.name = "Coefficient"
)

ggplot(df_coef_long, aes(x = Coefficient, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "Coefficient Distributions",
       x = "Coefficient", y = "Frequency")

# ---- PLOT coefficient differences ----
df_diff_long <- melt(
  data.frame(
    Diff_Coef_Maths  = results$Diff_Coef_Maths,
    Diff_Coef_Stoat  = results$Diff_Coef_Stoat,
    Diff_Coef_RVTest = results$Diff_Coef_RVTest
  ),
  variable.name = "Method",
  value.name = "CoefDiff"
)

ggplot(df_diff_long, aes(x = CoefDiff, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "Difference in Coefficients vs R lm",
       x = "C++ Coefficient - R Coefficient", y = "Frequency")

# ---- SIGNIFICANCE proportion ----
p_threshold <- 0.01
sig_props <- sapply(results[, c("P_R", "P_Maths", "P_Stoat", "P_RVTest")],
                    function(p) mean(p < p_threshold, na.rm = TRUE) * 100)

print(sig_props)


for (rep in 1:n_reps) {
  
  # Random sample sizes
  n <- sample(100:1000, 1)
  k <- 5 # total predictors incl. intercept [default : 2 path] + intercept
  number_last_element <- 1  # number of rows to modify at the end
  
  # Step 1: Generate matrix X (SNP data with sum=1 per row)
  X <- matrix(
    t(replicate(n, gen_row(k))),
    nrow = n,
    ncol = k
  )
  
  # Step 2: Add 2 new zero columns
  X <- cbind(X, matrix(0, nrow = n, ncol = 2))
  
  # Step 3: Modify the last 'number_last_element' rows
  if (number_last_element > 0) {
    start_row <- n - number_last_element + 1
    for (i in start_row:n) {
      X[i, ] <- 0                   # reset the whole row to 0
      X[i, (k+1):(k+2)] <- 0.5      # set last 2 columns to 0.5
    }
  }
  
  # Step 3: sanity check
  if(any(rowSums(X) != 1)) stop("Some rows do not sum to 1!")

  # ---- Store meta ----
  results$Rep[row_index] <- rep
  results$N[row_index]   <- n
  results$K[row_index]   <- k
  
  # ---- Generate Y ----
  beta_true <- c(2, rep(0, ncol(X)-1))
  intercept <- 1
  Y <- intercept + X %*% beta_true + rnorm(n, sd = 0.1)

  # ---- R lm ----
  df <- data.frame(Y = Y, X)
  fit_r <- lm(Y ~ ., data = df)
  coefs <- coef(summary(fit_r))
  
  results$P_R[row_index]     <- coefs[-1, "Pr(>|t|)"][1]
  results$Coef_R[row_index]  <- coefs[-1, "Estimate"][1]
  
  dup_cols <- duplicated(as.data.frame(t(X)))
  X_unique <- X[, !dup_cols, drop = FALSE]
  
  # Now proceed with maths version using X_unique instead of X
  X_maths <- X_unique[, -ncol(X_unique), drop = FALSE]
  
  # ---- C++ Maths ----
  res_maths <- cpp_linear_regression_maths(X_maths, Y)
  cpp_pvals_maths <- unlist(res_maths$p_values)
  results$P_Maths[row_index]     <- cpp_pvals_maths[2]
  results$Coef_Maths[row_index]  <- unlist(res_maths$coefficients)[2]
  results$Diff_Maths[row_index]  <- results$P_Maths[row_index] - results$P_R[row_index]
  
  # ---- C++ Stoat ----
  res_stoat <- cpp_linear_regression_stoat(X_maths, Y)
  cpp_pvals_stoat <- unlist(res_stoat$p_values)
  results$P_Stoat[row_index]     <- cpp_pvals_stoat[2]
  results$Coef_Stoat[row_index]  <- unlist(res_stoat$coefficients)[2]
  results$Diff_Stoat[row_index]  <- results$P_Stoat[row_index] - results$P_R[row_index]
  
  # ---- C++ RVTest ----
  res_rvtest <- cpp_linear_regression_rvtest(X_maths, Y)
  cpp_pvals_rvtest <- unlist(res_rvtest$p_values)
  results$P_RVTest[row_index]    <- cpp_pvals_rvtest[2]
  results$Coef_RVTest[row_index] <- unlist(res_rvtest$coefficients)[2]
  results$Diff_RVTest[row_index] <- results$P_RVTest[row_index] - results$P_R[row_index]
  
  # ---- Increment row index ----
  row_index <- row_index + 1
}

# ----------------------------------------- SIMULATION -----------------------------------------

# ---- Filtered p-value distributions: 1e-5 to 0.01 ----
df_pvals_long_small <- subset(df_pvals_long, PValue >= 1e-5 & PValue <= 0.01)

ggplot(df_pvals_long_small, aes(x = PValue, fill = Method)) +
  geom_histogram(alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "P-value Distributions by Method (1e-5 to 0.01)",
       x = "P-value", y = "Frequency") +
  xlim(1e-5, 0.01)

# ---- PLOT p-value distributions ----
df_pvals_long <- melt(
  data.frame(
    P_R      = results$P_R,
    P_Maths  = results$P_Maths,
    P_Stoat  = results$P_Stoat,
    P_RVTest = results$P_RVTest
  ),
  variable.name = "Method",
  value.name = "PValue"
)

ggplot(df_pvals_long, aes(x = PValue, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "P-value Distributions by Method",
       x = "P-value", y = "Frequency")

# ---- PLOT p-value differences ----
df_pdiff_long <- melt(
  data.frame(
    Diff_P_Maths  = results$Diff_Maths,
    Diff_P_Stoat  = results$Diff_Stoat,
    Diff_P_RVTest = results$Diff_RVTest
  ),
  variable.name = "Method",
  value.name = "PValueDiff"
)

ggplot(df_pdiff_long, aes(x = PValueDiff, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "Difference in P-values vs R lm",
       x = "C++ P-value - R P-value", y = "Frequency")

# ---- PLOT coefficient distributions ----
df_coef_long <- melt(
  data.frame(
    Coef_R      = results$Coef_R,
    Coef_Maths  = results$Coef_Maths,
    Coef_Stoat  = results$Coef_Stoat,
    Coef_RVTest = results$Coef_RVTest
  ),
  variable.name = "Method",
  value.name = "Coefficient"
)

ggplot(df_coef_long, aes(x = Coefficient, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "Coefficient Distributions",
       x = "Coefficient", y = "Frequency")

# ---- PLOT coefficient differences ----
df_diff_long <- melt(
  data.frame(
    Diff_Coef_Maths  = results$Diff_Coef_Maths,
    Diff_Coef_Stoat  = results$Diff_Coef_Stoat,
    Diff_Coef_RVTest = results$Diff_Coef_RVTest
  ),
  variable.name = "Method",
  value.name = "CoefDiff"
)

ggplot(df_diff_long, aes(x = CoefDiff, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "Difference in Coefficients vs R lm",
       x = "C++ Coefficient - R Coefficient", y = "Frequency")

# ---- SIGNIFICANCE proportion ----
p_threshold <- 0.01
sig_props <- sapply(results[, c("P_R", "P_Maths", "P_Stoat", "P_RVTest")],
                    function(p) mean(p < p_threshold, na.rm = TRUE) * 100)

print(sig_props)
