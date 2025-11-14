suppressPackageStartupMessages({
  library(car)
  library(dplyr)
  library(tibble)
})

# -------------------------------
# Dataset
# -------------------------------
y <- c(10.5, 13.0, 15.8, 19.7, 21.5)
x1 <- c(0.5, 0, 1, 0, 0)
x2 <- c(0, 0.5, 0, 1, 0.5)
cov1 <- c(1.1, 1.2, 0.8, 0.9, 1.0)

data <- tibble(
  y = y,
  x1 = x1,
  x2 = x2,
  cov1 = cov1
)

# -------------------------------
# R implementation (lm)
# -------------------------------
# Fit linear model including 2 predictors + covariate
fit_r <- lm(y ~ x1 + x2 + cov1, data = data)

# Extract global F-test from summary
summary_r <- summary(fit_r)
F_value <- as.numeric(summary_r$fstatistic["value"])
df1 <- as.numeric(summary_r$fstatistic["numdf"])
df2 <- as.numeric(summary_r$fstatistic["dendf"])
p_ftest_lm <- pf(F_value, df1, df2, lower.tail = FALSE)

# R² value
r2_lm <- summary_r$r.squared

# Linear hypothesis test for x1, x2, cov1
lh <- linearHypothesis(fit_r, c("x1=0", "x2=0"))
p_ftest_car <- lh[["Pr(>F)"]][2]

# -------------------------------
# Print results
# -------------------------------
cat("🔹 Global F-test (lm):\n")
cat("   p-value =", p_ftest_lm, "\n")
cat("   R² =", r2_lm, "\n\n")

cat("🔹 linearHypothesis (car):\n")
cat("   p-value =", p_ftest_car, "\n")
cat("   R² =", r2_lm, "\n")
