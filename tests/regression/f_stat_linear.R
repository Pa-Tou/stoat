# ---- Data ----
X <- matrix(c(
0,0,0,0,0,1,
1,1,1,0,0,2,
0,0,0,0,0,1,
0,1,1,0,0,1,
1,0,0,1,0,2,
0,0,0,0,1,2,
0,1,0,0,0,1,
0,0,0,0,0,1
), nrow=8, byrow=TRUE)

Y <- c(39.347,53.6,60.8,54.596,70.351,51.429,64.16,44.082)

# convert to dataframe
df <- data.frame(Y, X)
colnames(df) <- c("Y","X1","X2","X3","X4","X5","X6")

# ---- Models ----
full_model <- lm(Y ~ . - 1 , data=df)
# Reduced model: remove X1
reduced_model <- lm(Y ~ . - X1, data=df)

# check design matrices
full_X <- model.matrix(full_model)
reduced_X <- model.matrix(reduced_model)

cat("Full model design matrix:\n")
print(full_X)

cat("\nReduced model design matrix:\n")
print(reduced_X)

# check column rank
cat("\nRank full model:", qr(full_X)$rank, "\n")
cat("Rank reduced model:", qr(reduced_X)$rank, "\n")

# ---- Predictions ----
pred_full <- predict(full_model)
pred_reduced <- predict(reduced_model)

# ---- SSE ----
SSE_full <- sum(residuals(full_model)^2)
SSE_reduced <- sum(residuals(reduced_model)^2)

# ---- Degrees of freedom ----
df_full <- df.residual(full_model)
df_reduced <- df.residual(reduced_model)

df_numerator <- df_reduced - df_full
df_denominator <- df_full

# ---- Print everything ----
cat("----- DATA -----\n")
print(X)
print(Y)

cat("\n----- MODEL SUMMARY (FULL) -----\n")
print(summary(full_model))
print(summary(reduced_model))

cat("\n----- SSE -----\n")
cat("SSE_full:", SSE_full,"\n")
cat("SSE_reduced:", SSE_reduced,"\n")

cat("\n----- DF -----\n")
cat("df_full:", df_full,"\n")
cat("df_reduced:", df_reduced,"\n")
cat("df_numerator:", df_numerator,"\n")
cat("df_denominator:", df_denominator,"\n")

cat("\n----- F TEST PARTS -----\n")
cat("numerator:", numerator,"\n")
cat("denominator:", denominator,"\n")
cat("F_stat:", F_stat,"\n")

# [Running] Rscript "/home/mbagarre/Bureau/stoat/tempCodeRunnerFile.R"
# Full model design matrix:
#   (Intercept) X1 X2 X3 X4 X5 X6
# 1           1  0  0  0  0  0  1
# 2           1  1  1  1  0  0  2
# 3           1  0  0  0  0  0  1
# 4           1  0  1  1  0  0  1
# 5           1  1  0  0  1  0  2
# 6           1  0  0  0  0  1  2
# 7           1  0  1  0  0  0  1
# 8           1  0  0  0  0  0  1
# attr(,"assign")
# [1] 0 1 2 3 4 5 6

# Reduced model design matrix:
#   (Intercept) X2 X3 X4 X5 X6
# 1           1  0  0  0  0  1
# 2           1  1  1  0  0  2
# 3           1  0  0  0  0  1
# 4           1  1  1  0  0  1
# 5           1  0  0  1  0  2
# 6           1  0  0  0  1  2
# 7           1  1  0  0  0  1
# 8           1  0  0  0  0  1
# attr(,"assign")
# [1] 0 1 2 3 4 5

# Rank full model: 6 
# Rank reduced model: 6 
# ----- DATA -----
#      [,1] [,2] [,3] [,4] [,5] [,6]
# [1,]    0    0    0    0    0    1
# [2,]    1    1    1    0    0    2
# [3,]    0    0    0    0    0    1
# [4,]    0    1    1    0    0    1
# [5,]    1    0    0    1    0    2
# [6,]    0    0    0    0    1    2
# [7,]    0    1    0    0    0    1
# [8,]    0    0    0    0    0    1
# [1] 39.347 53.600 60.800 54.596 70.351 51.429 64.160 44.082

# ----- MODEL SUMMARY (FULL) -----

# Call:
# lm(formula = Y ~ ., data = df)

# Residuals:
#          1          2          3          4          5          6          7 
# -8.729e+00 -4.441e-16  1.272e+01  4.441e-16  8.882e-16  4.441e-16  2.220e-15 
#          8 
# -3.994e+00 

# Coefficients: (1 not defined because of singularities)
#             Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   48.076      6.507   7.388   0.0178 *
# X1            -0.996     15.939  -0.062   0.9559  
# X2            16.084     13.014   1.236   0.3420  
# X3            -9.564     15.939  -0.600   0.6094  
# X4            23.271     20.577   1.131   0.3755  
# X5             3.353     13.014   0.258   0.8208  
# X6                NA         NA      NA       NA  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 11.27 on 2 degrees of freedom
# Multiple R-squared:  0.6529,	Adjusted R-squared:  -0.2148 
# F-statistic: 0.7525 on 5 and 2 DF,  p-value: 0.6555


# ----- SSE -----
# SSE_full: 254.0477 
# SSE_reduced: 254.0477 

# ----- DF -----
# df_full: 2 
# df_reduced: 2 
# df_numerator: 0 
# df_denominator: 2