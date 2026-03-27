library(lmtest)
library(logistf)

# ----- Data -----

mat <- matrix(c(
1,0,1,2,
1,2,0,2,
1,0,2,2,
0,0,0,2,
2,0,0,1,
2,0,2,2,
0,0,0,2
), byrow=TRUE, ncol=4)

colnames(mat) <- paste0("A",0:3)
rownames(mat) <- paste0("S",1:7)

df <- as.data.frame(mat)

# binary phenotype
df$y <- as.logical(c(1,0,1,1,0,1,1))

# ----- Logistic regression -----

# full model
mod.f <- logistf(y~., data=df)

print("mod.f: ")
summary(mod.f)

# reduced model (no predictor)
mod.n <- logistf(y~1+A1+A2+A3, data=df)

print("mod.n: ")
summary(mod.n)

lr.chi2 = -2*(mod.n$loglik[1]-mod.f$loglik[1])

print(pchisq(lr.chi2, 1, lower.tail=FALSE))

# ----- Data -----

mat <- matrix(c(
2,0,0,1,2,
2,1,1,0,2,
1,0,2,0,2,
1,0,0,2,2,
0,0,0,0,1,
1,0,0,2,2,
2,0,0,0,2
), byrow=TRUE, ncol=5)

colnames(mat) <- paste0("A",0:4)
rownames(mat) <- paste0("S",1:7)

df <- as.data.frame(mat)

# binary phenotype
df$y <- as.logical(c(1,1,0,1,0,1,1))

# ----- Logistic regression -----

mod.f <- logistf(y~., data=df)

print("mod.f: ")
summary(mod.f)

# reduced model (no predictor)
mod.n <- logistf(y~1+A1+A2+A3+A4, data=df)

print("mod.n: ")
summary(mod.n)

lr.chi2 = -2*(mod.n$loglik[1]-mod.f$loglik[1])

print(pchisq(lr.chi2, 1, lower.tail=FALSE))

# [1] "mod.f: "
# logistf(formula = y ~ ., data = df)

# Model fitted by Penalized ML
# Coefficients:
#                   coef  se(coef) lower 0.95 upper 0.95      Chisq         p
# (Intercept) -5.0148450 10.009020 -34.362348  17.952131 0.24828108 0.6182880
# A0           0.2504847  2.366490  -5.348080   7.636251 0.01099036 0.9165068
# A1          -1.5823893  1.595444  -5.336683   1.684762 1.01083212 0.3147036
# A2          -0.4837771  1.969460  -4.500670   4.768237 0.05839643 0.8090486
# A3           3.4152634  5.129728  -8.203650  17.869491 0.43720481 0.5084747
#             method
# (Intercept)      2
# A0               2
# A1               2
# A2               2
# A3               2

# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None

# Likelihood ratio test=3.754772 on 4 df, p=0.4402098, n=7
# Wald test = 3.502074 on 4 df, p = 0.4775631[1] "mod.n: "
# logistf(formula = y ~ 1 + A1 + A2 + A3, data = df)

# Model fitted by Penalized ML
# Coefficients:
#                      coef se(coef) lower 0.95 upper 0.95        Chisq         p
# (Intercept) -3.988984e+00 3.638871 -14.369744  2.4869232 1.428232e+00 0.2320533
# A1          -1.445186e+00 1.144703  -4.875156  0.5284027 1.976933e+00 0.1597138
# A2           1.567930e-15 1.186611  -3.609160  3.6091604 6.217249e-15 0.9999999
# A3           2.890372e+00 2.289406  -1.056805  9.7503113 1.976933e+00 0.1597138
#             method
# (Intercept)      2
# A1               2
# A2               2
# A3               2

# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None

# Likelihood ratio test=4.180289 on 3 df, p=0.2426427, n=7
# Wald test = 3.656985 on 3 df, p = 0.3009652full 
#    1 
# [1] "mod.f: "
# logistf(formula = y ~ ., data = df)

# Model fitted by Penalized ML
# Coefficients:
#                      coef se(coef) lower 0.95 upper 0.95      Chisq         p
# (Intercept) -4.317488e+00 9.323804 -24.257496  18.061148 0.20841475 0.6480128
# A0          -5.108256e-01 3.966527  -8.446965  10.022757 0.01632581 0.8983290
# A1           1.354025e+00 2.569046  -4.480068   7.268358 0.27510588 0.5999271
# A2          -1.354025e+00 2.569046  -7.268358   4.480068 0.27510588 0.5999271
# A3          -1.003135e-14 2.309401  -5.552183   5.552183 0.00000000 1.0000000
# A4           3.218876e+00 8.884443 -18.715044  22.343784 0.12780869 0.7207150
#             method
# (Intercept)      2
# A0               2
# A1               2
# A2               2
# A3               2
# A4               2

# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None

# Likelihood ratio test=3.373169 on 5 df, p=0.6426595, n=7
# Wald test = 3.342317 on 5 df, p = 0.6473694[1] "mod.n: "
# logistf(formula = y ~ 1 + A1 + A2 + A3 + A4, data = df)

# Model fitted by Penalized ML
# Coefficients:
#                   coef se(coef) lower 0.95 upper 0.95      Chisq         p
# (Intercept) -3.5463215 3.687855 -13.992709  3.0361708 1.07791687 0.2991640
# A1           0.9733699 2.016615  -3.348015  6.3140374 0.24190457 0.6228339
# A2          -1.2238546 1.183245  -4.964895  0.8540869 1.25927815 0.2617878
# A3           0.2332924 1.243522  -3.738724  4.1102071 0.03201088 0.8580038
# A4           2.4477093 2.366490  -1.708174  9.9297902 1.25927815 0.2617878
#             method
# (Intercept)      2
# A1               2
# A2               2
# A3               2
# A4               2

# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None

# Likelihood ratio test=3.754772 on 4 df, p=0.4402098, n=7
# Wald test = 3.502074 on 4 df, p = 0.4775631full 
#    1