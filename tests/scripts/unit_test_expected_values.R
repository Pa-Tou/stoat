library(car)
library(logistf)
library(lmtest)

## Full dataset (same as C++)
X <- matrix(c(
     0, 0.5, 0, 10,  2,  5,
     1, 0,   0,  8,  1,  3,
     0, 1,   0,  9,  4,  2,
     1, 0,   0,  7,  3,  1,
     0, 0.5, 0, 11,  2,  4,
     1, 0,   0,  6,  1,  2,
     0, 1,   0, 10,  3,  5,
     1, 0,   0,  9,  2,  1,
     0, 0,   0,  8,  4,  3,
     0.5, 0, 0.5, 7,  2,  4
), byrow = TRUE, ncol = 6)

## Keep only last 3 columns as covariates
df <- data.frame(
    V5 = X[,4],
    V6 = X[,5],
    V7 = X[,6]
)

## Outcome (from C++)
y <- c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0)

df$y <- y

## Reduced logistic model (intercept is added automatically)
model_reduced <- glm(y ~ V5 + V6 + V7,
                     family = binomial,
                     data = df)

## Firth reduced model
fit_reduced <- logistf(y ~ V5 + V6 + V7, data = df)
fit_null <- logistf(y ~ 1, data = df)

## Likelihood ratio test (Firth)
lr_chi2 <- -2 * (fit_null$loglik[1] - fit_reduced$loglik[1])
pchisq(lr_chi2, df = 3, lower.tail = FALSE)

# p-value : 0.1667841
