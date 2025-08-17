# install.packages("RcppEigen", repos = "https://cloud.r-project.org")   # if not installed
# install.packages("RcppGSL", repos = "https://cloud.r-project.org")     # if not installed

library(Rcpp)
library(RcppEigen)
# library(RcppGSL)
library(ggplot2)
library(reshape2)

Sys.setenv("CXX14FLAGS"="-std=gnu++14")

# Load C++ implementations
sourceCpp("/Users/matisalias/Desktop/stoat/tests/regression/linear_regression_maths.cpp")
sourceCpp("/Users/matisalias/Desktop/stoat/tests/regression/linear_regression_stoat.cpp")
sourceCpp("/Users/matisalias/Desktop/stoat/tests/regression/linear_regression_rvtest.cpp")

# Simulation parameters
set.seed(123)
n_sims <- 1000
vals <- c(0, 0.5, 1)
p_threshold <- 0.01
intercept <- 1
sigma <- 1 # noise level

# Storage for simulation results
results <- data.frame(
  P_R = numeric(n_sims),
  P_Maths = numeric(n_sims),
  P_Stoat = numeric(n_sims),
  P_RVTest = numeric(n_sims),
  Diff_P_Maths = numeric(n_sims),
  Diff_P_Stoat = numeric(n_sims),
  Diff_P_RVTest = numeric(n_sims),
  N = integer(n_sims),
  K = integer(n_sims)
)

# Create storage for null results
results_null <- data.frame(
  P_R = numeric(n_sims),
  P_Maths = numeric(n_sims),
  P_Stoat = numeric(n_sims),
  P_RVTest = numeric(n_sims),
  Diff_P_Maths = numeric(n_sims),
  Diff_P_Stoat = numeric(n_sims),
  Diff_P_RVTest = numeric(n_sims)
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
  n <- 1000 # sample number
  k <- 3 # total predictors incl. intercept [default : 2 path] + intercept

  results$N[i] <- n
  results$K[i] <- k

  # Step 1a: Generate matrix X (SNP data with sum=1 per row)
  X <- matrix(
    t(replicate(n, gen_row(k-1))),
    nrow = n,
    ncol = k-1
  )

  # Step 1b: sanity check
  if(any(rowSums(X) != 1)) stop("Some rows do not sum to 1!")

  # ---- Step 2a: Significant regression ----
  beta_true <- rep(0, ncol(X))
  beta_true[1] <- 0.5  # assign signal to X1
  Y <- X %*% beta_true + rnorm(n, sd = sigma)
  Y <- as.vector(Y)

  # non-significant regression
  beta_null <- rep(0, ncol(X))
  beta_null[1] <- 0  # assign signal to X1
  Y_null <- X %*% beta_null + rnorm(n, sd = sigma)
  Y_null <- as.vector(Y_null)

  # ----------------------------------------- SIGNIFICATIVE -----------------------------------------

  # Step 3a: Linear regression in R
  df <- data.frame(Y = Y, X)
  fit_r <- lm(Y ~ ., data = df)
  r_summary <- summary(fit_r)
  coefs <- coef(r_summary)

  # Store the coefficient and p-value of the first predictor
  results$Coef_R[i] <- coef(fit_r)[2]
  results$P_R[i] <- coefs[-1, "Pr(>|t|)"][1] # Exclude intercept p-value

  # If p-value is NA, print summary and exit
  if (is.na(results$P_R[i])) {
    print(r_summary)  # Print full model summary
    stop("P-value is NA. Exiting.")
  }

  # exclude last column
  X_maths <- X[, -ncol(X), drop = FALSE]

  # Step 4: Linear regression using C++ maths implementation
  res_maths <- cpp_linear_regression_maths(X_maths, Y)
  cpp_coef_maths <- unlist(res_maths$coefficients)
  cpp_pvals_maths <- unlist(res_maths$p_values)

  results$P_Maths[i] <- cpp_pvals_maths[2]
  results$Coef_Maths[i] <- cpp_coef_maths[2]

  results$Diff_P_Maths[i] <- results$P_Maths[i] - results$P_R[i]
  results$Diff_Coef_Maths[i] <- results$Coef_Maths[i] - results$Coef_R[i]

  # Step 5: Linear regression using C++ stoat implementation
  res_stoat <- cpp_linear_regression_stoat(X_maths, Y)
  cpp_coef_stoat <- unlist(res_stoat$coefficients)
  cpp_pvals_stoat <- unlist(res_stoat$p_values)

  results$P_Stoat[i] <- cpp_pvals_stoat[2]
  results$Coef_Stoat[i] <- cpp_coef_stoat[2]

  results$Diff_P_Stoat[i]   <- results$P_Stoat[i] - results$P_R[i]
  results$Diff_Coef_Stoat[i] <- results$Coef_Stoat[i] - results$Coef_R[i]

  # Step 6: Linear regression using C++ rvtest implementation
  res_rvtest <- cpp_linear_regression_rvtest(X_maths, Y)
  cpp_coef_rvtest <- unlist(res_rvtest$coefficients)
  cpp_pvals_rvtest <- unlist(res_rvtest$p_values)

  results$P_RVTest[i] <- cpp_pvals_rvtest[2]
  results$Coef_RVTest[i] <- cpp_coef_rvtest[2]

  results$Diff_P_RVTest[i] <- results$P_RVTest[i] - results$P_R[i]
  results$Diff_Coef_RVTest[i] <- results$Coef_RVTest[i] - results$Coef_R[i]

  # ----------------------------------------- NO SIGNIFICATIVE -----------------------------------------

  df_null <- data.frame(Y = Y_null, X)
  fit_r <- lm(Y ~ ., data = df_null)
  r_summary <- summary(fit_r)

  results_null$P_R[i] <- coef(r_summary)[2, "Pr(>|t|)"]
  results_null$Coef_R[i] <- coef(fit_r)[2]

  # If p-value is NA, print summary and exit
  if (is.na(results_null$P_R[i])) {
    print(r_summary)  # Print full model summary
    stop("P-value is NA. Exiting.")
  }

  # C++ Maths
  res_maths <- cpp_linear_regression_maths(X[, -ncol(X), drop=FALSE], Y_null)
  results_null$P_Maths[i] <- res_maths$p_values[2]
  results_null$Coef_Maths[i] <- res_maths$coefficients[2]
  results_null$Diff_P_Maths[i] <- results_null$P_Maths[i] - results_null$P_R[i]

  # C++ Stoat
  res_stoat <- cpp_linear_regression_stoat(X[, -ncol(X), drop=FALSE], Y_null)
  results_null$P_Stoat[i] <- res_stoat$p_values[2]
  results_null$Coef_Stoat[i] <- res_stoat$coefficients[2]
  results_null$Diff_P_Stoat[i] <- results_null$P_Stoat[i] - results_null$P_R[i]

  # C++ RVTest
  res_rvtest <- cpp_linear_regression_rvtest(X[, -ncol(X), drop=FALSE], Y_null)
  results_null$P_RVTest[i] <- res_rvtest$p_values[2]
  results_null$Coef_RVTest[i] <- res_rvtest$coefficients[2]
  results_null$Diff_P_RVTest[i] <- results_null$P_RVTest[i] - results_null$P_R[i]

  # all other
  results_null$Diff_P_Maths[i] <- results_null$P_Maths[i] - results_null$P_R[i]
  results_null$Diff_P_Stoat[i] <- results_null$P_Stoat[i] - results_null$P_R[i]
  results_null$Diff_P_RVTest[i]  <- results_null$P_RVTest[i] - results_null$P_R[i]

  results_null$Diff_Coef_Maths[i] <- results_null$Coef_Maths[i] - results_null$Coef_R[i]
  results_null$Diff_Coef_Stoat[i] <- results_null$Coef_Stoat[i] - results_null$Coef_R[i]
  results_null$Diff_Coef_RVTest[i] <- results_null$Coef_RVTest[i] - results_null$Coef_R[i]
}

# -------------------------------------------------------- PLOT SIGNIFICANCE ----------------------------------------------------

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
  labs(title = "P-value Distributions by Method [all significative]",
       x = "P-value", y = "Frequency")

# ---- Filtered p-value distributions: 0.01 to 1e-5 ----
df_pvals_long_small <- subset(df_pvals_long, PValue >= 1e-8 & PValue <= 0.01)

ggplot(df_pvals_long_small, aes(x = PValue, fill = Method)) +
  geom_histogram(alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "P-value Distributions by Method (0.01 to 1e-8) [all significative]",
       x = "P-value", y = "Frequency") +
  xlim(1e-8, 0.01)

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
  labs(title = "Difference in P-values vs R lm [all significative]",
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
  labs(title = "Coefficient Distributions [all significative]",
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
  labs(title = "Difference in Coefficients vs R lm [all significative]",
       x = "C++ Coefficient - R Coefficient", y = "Frequency")

# ---- SIGNIFICANCE proportion ----
p_threshold <- 0.01
sig_props <- sapply(results[, c("P_R", "P_Maths", "P_Stoat", "P_RVTest")],
                    function(p) mean(p < p_threshold, na.rm = TRUE) * 100)

print(sig_props)

# ---- Difference proportion ----
cat("Difference Pvalue per methods [min/max/means] : \n")

diff_stats <- data.frame(
  Method = c("Maths", "Stoat", "RVTest"),
  Min = c(
    min(results$Diff_P_Maths, na.rm = TRUE),
    min(results$Diff_P_Stoat, na.rm = TRUE),
    min(results$Diff_P_RVTest, na.rm = TRUE)
  ),
  Mean = c(
    mean(results$Diff_P_Maths, na.rm = TRUE),
    mean(results$Diff_P_Stoat, na.rm = TRUE),
    mean(results$Diff_P_RVTest, na.rm = TRUE)
  ),
  Max = c(
    max(results$Diff_P_Maths, na.rm = TRUE),
    max(results$Diff_P_Stoat, na.rm = TRUE),
    max(results$Diff_P_RVTest, na.rm = TRUE)
  )
)

print(diff_stats)

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
  labs(title = "P-value Distributions by Method [all NO significative]",
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
  labs(title = "Difference in P-values vs R lm [all NO significative]",
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
  labs(title = "Coefficient Distributions [all NO significative]",
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
  labs(title = "Difference in Coefficients vs R lm [all NO significative]",
       x = "C++ Coefficient - R Coefficient", y = "Frequency")

# ---- SIGNIFICANCE proportion for non-significant simulations ----
p_threshold <- 0.01

sig_props_null <- sapply(
  results_null[, c("P_R", "P_Maths", "P_Stoat", "P_RVTest")],
  function(p) mean(p < p_threshold, na.rm = TRUE) * 100
)

print(sig_props_null)

# ---- Difference proportion ----
cat("Difference Pvalue per methods [min/max/means] : \n")

diff_stats_null <- data.frame(
  Method = c("Maths", "Stoat", "RVTest"),
  Min = c(
    min(results_null$Diff_P_Maths, na.rm = TRUE),
    min(results_null$Diff_P_Stoat, na.rm = TRUE),
    min(results_null$Diff_P_RVTest, na.rm = TRUE)
  ),
  Mean = c(
    mean(results_null$Diff_P_Maths, na.rm = TRUE),
    mean(results_null$Diff_P_Stoat, na.rm = TRUE),
    mean(results_null$Diff_P_RVTest, na.rm = TRUE)
  ),
  Max = c(
    max(results_null$Diff_P_Maths, na.rm = TRUE),
    max(results_null$Diff_P_Stoat, na.rm = TRUE),
    max(results_null$Diff_P_RVTest, na.rm = TRUE)
  )
)

print(diff_stats_null)

# ----------------------------------------- SIMULATION COLLINEARITY -----------------------------------------

# ----------------------------------------- SIGNIFICATIF -----------------------------------------

results <- data.frame(
  Rep   = integer(n_sims),
  N     = integer(n_sims),
  K     = integer(n_sims),
  P_R = numeric(n_sims),
  P_Maths = numeric(n_sims),
  P_Stoat = numeric(n_sims),
  P_RVTest = numeric(n_sims),
  Coef_R = numeric(n_sims),
  Coef_Maths = numeric(n_sims),
  Coef_Stoat = numeric(n_sims),
  Coef_RVTest = numeric(n_sims),
  Diff_P_Maths = numeric(n_sims),
  Diff_P_Stoat = numeric(n_sims),
  Diff_P_RVTest = numeric(n_sims),
  Diff_Coef_Maths  <- numeric(n_sims),
  Diff_Coef_Stoat  <- numeric(n_sims),
  Diff_Coef_RVTest <- numeric(n_sims)
)

for (rep in 1:n_sims) {

  # Random sample sizes
  n <- 1000
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
  results$Rep[rep] <- rep
  results$N[rep]   <- n
  results$K[rep]   <- k

  # ---- Generate Y ----
  beta_true <- rep(0, ncol(X))
  beta_true[1] <- k  # assign signal to X1
  Y <- X %*% beta_true + rnorm(n, sd = sigma)
  Y <- as.vector(Y)

  # ---- R lm ----
  df <- data.frame(Y = Y, X)
  fit_r <- lm(Y ~ ., data = df)
  coefs <- coef(summary(fit_r))

  results$P_R[rep]     <- coefs[2, "Pr(>|t|)"]
  results$Coef_R[rep]  <- coef(fit_r)[2]

  dup_cols <- duplicated(as.data.frame(t(X)))
  X_unique <- X[, !dup_cols, drop = FALSE]

  # Now proceed with maths version using X_unique instead of X
  X_maths <- X_unique[, -ncol(X_unique), drop = FALSE]

  # ---- C++ Maths ----
  res_maths <- cpp_linear_regression_maths(X_maths, Y)
  cpp_pvals_maths <- unlist(res_maths$p_values)
  results$P_Maths[rep]     <- cpp_pvals_maths[2]
  results$Coef_Maths[rep]  <- unlist(res_maths$coefficients)[2]
  results$Diff_P_Maths[rep]  <- results$P_Maths[rep] - results$P_R[rep]

  # ---- C++ Stoat ----
  res_stoat <- cpp_linear_regression_stoat(X_maths, Y)
  cpp_pvals_stoat <- unlist(res_stoat$p_values)
  results$P_Stoat[rep]     <- cpp_pvals_stoat[2]
  results$Coef_Stoat[rep]  <- unlist(res_stoat$coefficients)[2]
  results$Diff_P_Stoat[rep]  <- results$P_Stoat[rep] - results$P_R[rep]

  # ---- C++ RVTest ----
  res_rvtest <- cpp_linear_regression_rvtest(X_maths, Y)
  cpp_pvals_rvtest <- unlist(res_rvtest$p_values)
  results$P_RVTest[rep]    <- cpp_pvals_rvtest[2]
  results$Coef_RVTest[rep] <- unlist(res_rvtest$coefficients)[2]
  results$Diff_P_RVTest[rep] <- results$P_RVTest[rep] - results$P_R[rep]

}

# ----------------------------------------- SIMULATION -----------------------------------------

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

# ---- Filtered p-value distributions: 1e-5 to 0.01 ----
df_pvals_long_small <- subset(df_pvals_long, PValue >= 1e-8 & PValue <= 0.01)

ggplot(df_pvals_long_small, aes(x = PValue, fill = Method)) +
  geom_histogram(alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "P-value Distributions by Method (1e-8 to 0.01) [collinearity significative]",
       x = "P-value", y = "Frequency") +
  xlim(1e-8, 0.01)

ggplot(df_pvals_long, aes(x = PValue, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "P-value Distributions by Method [collinearity significative]",
       x = "P-value", y = "Frequency")

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
  labs(title = "Difference in P-values vs R lm [collinearity significative]",
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
  labs(title = "Coefficient Distributions [collinearity significative]",
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
  labs(title = "Difference in Coefficients vs R lm [collinearity significative]",
       x = "C++ Coefficient - R Coefficient", y = "Frequency")

# ---- SIGNIFICANCE proportion ----
p_threshold <- 0.01
sig_props <- sapply(results[, c("P_R", "P_Maths", "P_Stoat", "P_RVTest")],
                    function(p) mean(p < p_threshold, na.rm = TRUE) * 100)

print(sig_props)

cat("Difference Pvalue per methods [min/max/means] : \n")

diff_stats <- data.frame(
  Method = c("Maths", "Stoat", "RVTest"),
  Min = c(
    min(results$Diff_P_Maths, na.rm = TRUE),
    min(results$Diff_P_Stoat, na.rm = TRUE),
    min(results$Diff_P_RVTest, na.rm = TRUE)
  ),
  Mean = c(
    mean(results$Diff_P_Maths, na.rm = TRUE),
    mean(results$Diff_P_Stoat, na.rm = TRUE),
    mean(results$Diff_P_RVTest, na.rm = TRUE)
  ),
  Max = c(
    max(results$Diff_P_Maths, na.rm = TRUE),
    max(results$Diff_P_Stoat, na.rm = TRUE),
    max(results$Diff_P_RVTest, na.rm = TRUE)
  )
)

print(diff_stats)

# ----------------------------------------- NO SIGNIFICATIF -----------------------------------------

results_null <- data.frame(
  Rep   = integer(n_sims),
  N     = integer(n_sims),
  K     = integer(n_sims),
  P_R = numeric(n_sims),
  P_Maths = numeric(n_sims),
  P_Stoat = numeric(n_sims),
  P_RVTest = numeric(n_sims),
  Coef_R = numeric(n_sims),
  Coef_Maths = numeric(n_sims),
  Coef_Stoat = numeric(n_sims),
  Coef_RVTest = numeric(n_sims),
  Diff_P_Maths = numeric(n_sims),
  Diff_P_Stoat = numeric(n_sims),
  Diff_P_RVTest = numeric(n_sims),
  Diff_Coef_Maths  <- numeric(n_sims),
  Diff_Coef_Stoat  <- numeric(n_sims),
  Diff_Coef_RVTest <- numeric(n_sims)
)

for (rep in 1:n_sims) {

  # Random sample sizes
  n <- 1000
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
  results_null$Rep[rep] <- rep
  results_null$N[rep]   <- n
  results_null$K[rep]   <- k

  # ---- Generate Y ----
  beta_null <- rep(0, ncol(X))
  beta_null[1] <- 0  # assign signal to X1
  Y_null <- X %*% beta_null + rnorm(n, sd = sigma)
  Y_null <- as.vector(Y_null)

  # ---- R lm ----
  df_null <- data.frame(Y = Y_null, X)
  fit_r <- lm(Y ~ ., data = df_null)
  coefs <- coef(summary(fit_r))

  results_null$P_R[rep]     <- coefs[2, "Pr(>|t|)"]
  results_null$Coef_R[rep]  <- coef(fit_r)[2]

  dup_cols <- duplicated(as.data.frame(t(X)))
  X_unique <- X[, !dup_cols, drop = FALSE]

  # Now proceed with maths version using X_unique instead of X
  X_maths <- X_unique[, -ncol(X_unique), drop = FALSE]

  # ---- C++ Maths ----
  res_maths <- cpp_linear_regression_maths(X_maths, Y_null)
  cpp_pvals_maths <- unlist(res_maths$p_values)
  results_null$P_Maths[rep]     <- cpp_pvals_maths[2]
  results_null$Coef_Maths[rep]  <- unlist(res_maths$coefficients)[2]
  results_null$Diff_P_Maths[rep]  <- results_null$P_Maths[rep] - results_null$P_R[rep]

  # ---- C++ Stoat ----
  res_stoat <- cpp_linear_regression_stoat(X_maths, Y_null)
  cpp_pvals_stoat <- unlist(res_stoat$p_values)
  results_null$P_Stoat[rep]     <- cpp_pvals_stoat[2]
  results_null$Coef_Stoat[rep]  <- unlist(res_stoat$coefficients)[2]
  results_null$Diff_P_Stoat[rep]  <- results_null$P_Stoat[rep] - results_null$P_R[rep]

  # ---- C++ RVTest ----
  res_rvtest <- cpp_linear_regression_rvtest(X_maths, Y_null)
  cpp_pvals_rvtest <- unlist(res_rvtest$p_values)
  results_null$P_RVTest[rep]    <- cpp_pvals_rvtest[2]
  results_null$Coef_RVTest[rep] <- unlist(res_rvtest$coefficients)[2]
  results_null$Diff_P_RVTest[rep] <- results_null$P_RVTest[rep] - results_null$P_R[rep]

}

# ----------------------------------------- SIMULATION -----------------------------------------

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

# ---- PLOT p-value distributions ----
df_pvals_long <- melt(
  data.frame(
    P_R      = results_null$P_R,
    P_Maths  = results_null$P_Maths,
    P_Stoat  = results_null$P_Stoat,
    P_RVTest = results_null$P_RVTest
  ),
  variable.name = "Method",
  value.name = "PValue"
)

ggplot(df_pvals_long, aes(x = PValue, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "P-value Distributions by Method [collinearity NO significative]",
       x = "P-value", y = "Frequency")

# ---- PLOT p-value differences ----
df_pdiff_long <- melt(
  data.frame(
    Diff_P_Maths  = results_null$Diff_P_Maths,
    Diff_P_Stoat  = results_null$Diff_P_Stoat,
    Diff_P_RVTest = results_null$Diff_P_RVTest
  ),
  variable.name = "Method",
  value.name = "PValueDiff"
)

ggplot(df_pdiff_long, aes(x = PValueDiff, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "Difference in P-values vs R lm [collinearity NO significative]",
       x = "C++ P-value - R P-value", y = "Frequency")

# ---- PLOT coefficient distributions ----
df_coef_long <- melt(
  data.frame(
    Coef_R      = results_null$Coef_R,
    Coef_Maths  = results_null$Coef_Maths,
    Coef_Stoat  = results_null$Coef_Stoat,
    Coef_RVTest = results_null$Coef_RVTest
  ),
  variable.name = "Method",
  value.name = "Coefficient"
)

ggplot(df_coef_long, aes(x = Coefficient, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "Coefficient Distributions [collinearity NO significative]",
       x = "Coefficient", y = "Frequency")

# ---- PLOT coefficient differences ----
df_diff_long <- melt(
  data.frame(
    Diff_Coef_Maths  = results_null$Diff_Coef_Maths,
    Diff_Coef_Stoat  = results_null$Diff_Coef_Stoat,
    Diff_Coef_RVTest = results_null$Diff_Coef_RVTest
  ),
  variable.name = "Method",
  value.name = "CoefDiff"
)

ggplot(df_diff_long, aes(x = CoefDiff, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "dodge") +
  theme_bw() +
  labs(title = "Difference in Coefficients vs R lm [collinearity NO significative]",
       x = "C++ Coefficient - R Coefficient", y = "Frequency")

# ---- SIGNIFICANCE proportion ----
p_threshold <- 0.01
sig_props <- sapply(results_null[, c("P_R", "P_Maths", "P_Stoat", "P_RVTest")],
                    function(p) mean(p < p_threshold, na.rm = TRUE) * 100)

print(sig_props)

cat("Difference Pvalue per methods [min/max/means] : \n")

test_PR <- data.frame(
  Method = c("PR"),
  Min = min(results_null$P_R, na.rm = TRUE),
  Mean = mean(results_null$P_R, na.rm = TRUE),
  Max = max(results_null$P_R, na.rm = TRUE)
)

print(test_PR)

diff_stats_null <- data.frame(
  Method = c("Maths", "Stoat", "RVTest"),
  Min = c(
    min(results_null$Diff_P_Maths, na.rm = TRUE),
    min(results_null$Diff_P_Stoat, na.rm = TRUE),
    min(results_null$Diff_P_RVTest, na.rm = TRUE)
  ),
  Mean = c(
    mean(results_null$Diff_P_Maths, na.rm = TRUE),
    mean(results_null$Diff_P_Stoat, na.rm = TRUE),
    mean(results_null$Diff_P_RVTest, na.rm = TRUE)
  ),
  Max = c(
    max(results_null$Diff_P_Maths, na.rm = TRUE),
    max(results_null$Diff_P_Stoat, na.rm = TRUE),
    max(results_null$Diff_P_RVTest, na.rm = TRUE)
  )
)

print(diff_stats_null)