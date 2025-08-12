# install.packages("RcppEigen")   # if not installed
# install.packages("RcppGSL")     # if not installed

library(Rcpp)
library(RcppEigen)
library(RcppGSL)
library(ggplot2)
library(reshape2)
library(glmnet)

# Load C++ implementations
sourceCpp("linear_regression_maths.cpp")
sourceCpp("linear_regression_stoat.cpp")
sourceCpp("linear_regression_rvtest.cpp")

# Simulation parameters
set.seed(123)
n_sims <- 1000  # try small first
vals <- c(0, 0.5, 1)
p_threshold <- 0.01

# Generate one constrained SNP row (sum of SNPs = 1)
gen_row <- function(num_snps) {
  repeat {
    row <- sample(vals, num_snps, replace = TRUE)
    if (sum(row) == 1) return(row)
  }
}

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
for (i in seq_len(n_sims)) {
  # Random sample sizes
  n <- sample(100:1000, 1)
  k <- sample(3:6, 1) # total predictors incl. intercept [default : 2 path min] + intercept
  
  results$N[i] <- n
  results$K[i] <- k
  
  # Step 1: Generate matrix X (SNP data with sum=1 per row)
  snp_matrix <- matrix(
    t(replicate(n, gen_row(k - 1))),
    nrow = n,
    ncol = k - 1
  )

  X <- cbind(1, snp_matrix) # Add intercept column
  
  # Step 2: Generate phenotype vector Y
  beta_true <- runif(k - 1, -2, 2) # true effect sizes 
  intercept <- runif(1, 3, 6)
  Y <- intercept + snp_matrix %*% beta_true + rnorm(n, sd = 0.5)
  Y <- as.numeric(Y)
  
  # Step 3a: Linear regression in R
  df <- data.frame(Y = Y, snp_matrix)
  fit_r <- lm(Y ~ ., data = df)
  r_summary <- summary(fit_r)
  coefs <- coef(r_summary)
  p_r <- coefs[-1, "Pr(>|t|)"]  # exclude intercept p-value
  
  # For simplicity, take the p-value of the first effect (assuming at least one)
  results$P_R[i] <- p_r[2]
  r_coef <- coef(fit_r)[2]
  
  # Prepare X_maths: exclude last column if more than 1 SNP column to match your cpp funcs
  if (ncol(snp_matrix) > 1) {
    X_maths <- snp_matrix[, -ncol(snp_matrix), drop = FALSE]
  } else {
    X_maths <- snp_matrix
  }

  # Step 3b: Ridge regression in R
  # glmnet requires matrix X and vector y, without intercept column
  # X_ridge <- as.matrix(snp_matrix)  # Keep all SNPs for ridge
  # y_ridge <- Y

  # ridge_fit <- cv.glmnet(X_ridge, y_ridge, alpha = 0)  # alpha=0 → ridge regression
  # best_lambda <- ridge_fit$lambda.min

  # # Extract coefficients at optimal lambda
  # ridge_coefs <- coef(ridge_fit, s = "lambda.min")

  # # Save the coefficient of the first SNP predictor
  # if (nrow(ridge_coefs) > 1) {
  #   results$Coef_Ridge[i] <- ridge_coefs[2, 1]  # SNP1 effect
  # } else {
  #   results$Coef_Ridge[i] <- NA
  # }

  # Step 4: Linear regression using C++ maths implementation
  res_maths <- cpp_linear_regression_maths(X_maths, Y)
  cpp_coef_maths <- unlist(res_maths$coefficients)
  cpp_pvals_maths <- unlist(res_maths$p_values)
  results$P_Maths[i] <- cpp_pvals_maths[2]
  results$Diff_Maths[i] <- cpp_coef_maths[2] - r_coef

  # Step 5: Linear regression using C++ stoat implementation
  res_stoat <- cpp_linear_regression_stoat(X_maths, Y)
  cpp_coef_stoat <- unlist(res_stoat$coefficients)
  cpp_pvals_stoat <- unlist(res_stoat$p_values)
  results$P_Stoat[i] <- cpp_pvals_stoat[2]
  results$Diff_Stoat[i] <- cpp_coef_stoat[2] - r_coef

  # Step 6: Linear regression using C++ rvtest implementation
  res_rvtest <- cpp_linear_regression_rvtest(X_maths, Y)
  cpp_coef_rvtest <- unlist(res_rvtest$coefficients)
  cpp_pvals_rvtest <- unlist(res_rvtest$p_values)
  results$P_RVTest[i] <- cpp_pvals_rvtest[2]
  results$Diff_RVTest[i] <- cpp_coef_rvtest[2] - r_coef

  # cat("P-value:\n")
  # cat("R lm: ", p_r[2], "\n")
  # cat("C++ Maths: ", cpp_pvals_maths[2], "\n")
  # cat("C++ Stoat: ", cpp_pvals_stoat[2], "\n")
  # cat("C++ RVTest: ", cpp_pvals_rvtest[2], "\n\n")

  # cat("Coefficients:\n")
  # cat("R lm coef: ", r_coef, "\n")
  # cat("C++ Maths coef (first): ", cpp_coef_maths[2], "\n")
  # cat("C++ Stoat coef (first): ", cpp_coef_stoat[2], "\n")
  # cat("C++ RVTest coef (first): ", cpp_coef_rvtest[2], "\n")
}

# ---- PLOT p-value distributions ----
library(reshape2)
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

results$P_R

# Find the rows with non-finite p-values
bad_pvals <- subset(df_pvals_long, !is.finite(PValue))
print(table(bad_pvals$Method))
ggplot(df_pvals_long, aes(x = PValue, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
  theme_bw() +
  labs(title = "P-value Distributions", x = "P-value", y = "Frequency")

# ---- PLOT coefficient differences ----
df_diff_long <- melt(
  data.frame(
    Diff_Maths  = results$Diff_Maths,
    Diff_Stoat  = results$Diff_Stoat,
    Diff_RVTest = results$Diff_RVTest
  ),
  variable.name = "Method",
  value.name = "CoefDiff"
)

ggplot(df_diff_long, aes(x = CoefDiff, fill = Method)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
  theme_bw() +
  labs(title = "Difference in Coefficients vs R lm",
       x = "C++ Coefficient - R Coefficient", y = "Frequency")

# ---- SIGNIFICANCE proportion ----
p_threshold <- 0.05
sig_props <- sapply(results[, c("P_R", "P_Maths", "P_Stoat", "P_RVTest")],
                    function(p) mean(p < p_threshold, na.rm = TRUE) * 100)

print(sig_props)
