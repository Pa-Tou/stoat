library(Rcpp)
library(ggplot2)

# Load C++ implementation
sourceCpp("linear_regression_maths.cpp")
sourceCpp("linear_regression_stoat.cpp")
sourceCpp("linear_regression_rvtest.cpp")

# Step 1: Simulate data
set.seed(123)
n <- 50
k <- 3
vals <- c(0, 0.5, 1)

gen_row <- function(k) {
  repeat {
    row <- sample(vals, k, replace = TRUE)
    if (sum(row) == 1) return(row)
  }
}
X <- t(replicate(n, gen_row(k)))

# Step 2: Quantitative phenotype
beta_true <- c(2, -1, 0.5)
intercept <- 5
Y <- intercept + X %*% beta_true + rnorm(n, sd = 0.5)
Y <- as.numeric(Y)

# Step 3: R lm()
fit_r <- lm(Y ~ X)
r_coef <- coef(fit_r)

# Step 4: C++ regression
cpp_coef_maths <- unlist(cpp_linear_regression_maths(X, Y)$coefficients)
names(cpp_coef_maths) <- c("(Intercept)", paste0("X", 1:k))

cpp_coef_stoat <- unlist(cpp_linear_regression_stoat(X, Y)$coefficients)
names(cpp_coef_stoat) <- c("(Intercept)", paste0("X", 1:k))

cpp_coef_rvtest <- unlist(cpp_linear_regression_rvtest(X, Y)$coefficients)
names(cpp_coef_rvtest) <- c("(Intercept)", paste0("X", 1:k))

# Step 5: Compare
df_compare <- data.frame(
  Term = names(r_coef),
  R_Coefficient = r_coef,
  CPP_Coefficient = cpp_coef
)

print(df_compare)

# Step 6: Plot
ggplot(df_compare, aes(x = R_Coefficient, y = CPP_Coefficient, label = Term)) +
  geom_point(color = "blue", size = 3) +
  geom_text(vjust = -1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "R vs C++ Linear Regression Coefficients",
       x = "R Coefficients",
       y = "C++ Coefficients")
