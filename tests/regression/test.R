# Load required package
library(car)

# Simulate data
set.seed(123)
n <- 100

# Covariates (independent variables)
x1 <- rnorm(n, mean = 5, sd = 2)
x2 <- rnorm(n, mean = 10, sd = 3)
x3 <- rnorm(n, mean = 20, sd = 4)
x4 <- rnorm(n, mean = 0, sd = 1)

# Control variable (a covariate)
age <- rnorm(n, mean = 40, sd = 10)

# Dependent variable influenced by x1, x2, and covariate age
y <- 3 + 1.5*x1 - 1.6*x2 + rnorm(n, sd = 5)

# Combine into a data frame
data <- data.frame(y, x1, x2, x3, x4)
data_2 <- data.frame(y, x1, x2, x3, x4, age)

# Fit linear regression model including the covariate
fit_r <- lm(y ~ x1 + x2 + x3 + x4, data = data)

# Show results
print(summary(fit_r))

# Test the joint hypothesis that x1, x2, x3, and x4 have no effect
lh.o <- linearHypothesis(fit_r, c("x1 = 0", "x2 = 0", "x3 = 0", "x4 = 0"))

# Show results
print(lh.o)

# Fit linear regression model including the covariate
fit_r <- lm(y ~ x1 + x2 + x3 + x4 + age, data = data_2)

# Show results
print(summary(fit_r))

# Test the joint hypothesis that x1, x2, x3, and x4 have no effect
lh.o <- linearHypothesis(fit_r, c("x1 = 0", "x2 = 0", "x3 = 0", "x4 = 0"))

# Show results
print(lh.o)
