suppressPackageStartupMessages({
library(Rcpp)
library(RcppEigen)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(parallel)
library(car)
})

Sys.setenv("CXX14FLAGS"="-std=gnu++14")

sourceCpp("/home/mbagarre/Bureau/stoat/tests/regression/linear_regression_stoat.cpp")

# -------------------------------
# Single highly significant test
# -------------------------------

# Simulate one dataset with a strong effect
set.seed(42)
N <- 500

data <- tibble(y = rnorm(N)) %>%
  mutate(
    # x1 strongly correlated with y → very significant
    x1 = rnorm(N, y, sd = 2),
    # other predictors random noise
    x2 = rnorm(N),
    x3 = rnorm(N),
    x4 = rnorm(N)
  )

# Prepare matrices for Stoat
X <- as.matrix(data[, c("x1", "x2", "x3", "x4")])
Y <- data$y
Xcovar <- matrix(nrow = nrow(data), ncol = 0)  # no covariates

# -------------------------------
# R implementation
# -------------------------------
fit_r <- lm(y ~ x1 + x2 + x3 + x4, data = data)

# Global F-test (from summary)
summary_r <- summary(fit_r)
F_value <- as.numeric(summary_r$fstatistic["value"])
df1 <- as.numeric(summary_r$fstatistic["numdf"])
df2 <- as.numeric(summary_r$fstatistic["dendf"])
p_ftest_lm <- pf(F_value, df1, df2, lower.tail = FALSE)

# Linear hypothesis test (x1..x4)
lh.o <- car::linearHypothesis(fit_r, c("x1=0", "x2=0", "x3=0", "x4=0"))
p_ftest_car <- lh.o[["Pr(>F)"]][2]

# -------------------------------
# Stoat implementation (C++)
# -------------------------------
res_stoat <- cpp_linear_regression_stoat(X, Xcovar, Y)
p_stoat <- res_stoat$p_value

# -------------------------------
# Compare results
# -------------------------------
cat("🔹 Global F-test (lm):", p_ftest_lm, "\n")
cat("🔹 linearHypothesis (car):", p_ftest_car, "\n")
cat("🔹 Stoat (C++):", p_stoat, "\n")

# Optional: show -log10 comparison
tibble(
  Method = c("lm", "car", "Stoat"),
  `-log10(p-value)` = -log10(c(p_ftest_lm, p_ftest_car, p_stoat))
)

