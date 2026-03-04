## We used this test to set the expected values in the unit tests of the statistical tests. 
## Basically, we used R functions and packages to run the Chi2 or regression tests (linear, logistic, and Firth logistic) and extracted the p-value that we should expect.

library(car) ## for the F-test excluding covariates
library(logistf) ## to test Firth regression as an more robust logistic regression
library(lmtest)

## Logistic Regression 2 paths - Perfect Relationship
mat = matrix(c(0, 1,
               1, 0,
               1, 0,
               1, 0,
               0, 1, 
               0, 1,
               1, 0,
               1, 0,
               1, 0,
               0, 1), byrow=TRUE, ncol=2)
colnames(mat) = c('A0', 'A1')
rownames(mat) = paste0("S", 1:10)
df = as.data.frame(mat)
df$y = as.logical(c(1, 0, 0, 0, 1, 1, 0, 0, 0, 1))

((model = glm(y~A1, family="binomial", data=df)))
summary(model)$coef

summary(logistf(y~A1, data=df))
mod.f = logistf(y~A1, data=df)
mod.n = logistf(y~1, data=df)
lr.chi2 = -2*(mod.n$loglik[1]-mod.f$loglik[1])
pchisq(lr.chi2, 1, lower.tail=FALSE)

## Logistic Regression 2 paths - Moderate
mat = matrix(c(0, 1,
               1, 0,
               1, 0,
               1, 0,
               0, 1), byrow=TRUE, ncol=2)
colnames(mat) = c('A0', 'A1')
rownames(mat) = paste0("S", 1:5)
df = as.data.frame(mat)
df$y = as.logical(c(0, 1, 0, 0, 1))

mod.f = glm(y~A1, family="binomial", data=df)
mod.n = glm(y~1, family="binomial", data=df)

lrtest(mod.f, mod.n)

summary(logistf(y~A1, data=df))

## Logistic Regression 3 paths - Moderate
mat = matrix(c(1, 0, 1,
               1, 0, 1,
               1, 1, 0,
               0, 1, 1,
               0, 1, 1), byrow=TRUE, ncol=3)
colnames(mat) = c('A0', 'A1', 'A2')
rownames(mat) = paste0("S", 1:5)
df = as.data.frame(mat)
df$y = as.logical(c(1, 1, 1, 0, 0))

summary(glm(y~A1+A2, family="binomial", data=df))

summary(logistf(y~A1+A2, data=df))

mod.f = logistf(y~A1+A2, data=df)
mod.n = logistf(y~1, data=df)
lr.chi2 = -2*(mod.n$loglik[1]-mod.f$loglik[1])
pchisq(lr.chi2, 2, lower.tail=FALSE)


## Linear Regression 2 paths - good linear relationship
mat = matrix(c(1, 0,
               1, 0,
               1, 0,
               0, 1,
               0, 1), byrow=TRUE, ncol=2)
colnames(mat) = c('A0', 'A1')
rownames(mat) = paste0("S", 1:5)
df = as.data.frame(mat)
df$y = c(11, 10.1, 5.2, -0.3, 2)

summary(lm(y~A1, data=df))

## Linear Regression 3 paths - Moderate
mat = matrix(c(1, 0, 1,
               1, 0, 1,
               1, 1, 0,
               0, 1, 1,
               0, 1, 1), byrow=TRUE, ncol=3)
colnames(mat) = c('A0', 'A1', 'A2')
rownames(mat) = paste0("S", 1:5)
df = as.data.frame(mat)
df$y = c(11, 10.1, 5.2, -0.3, 2)

summary(lm(y~A1+A2, data=df))

## Linear Regression pseudo inversion 3 paths - Moderate
mat = matrix(c(1, 0, 0,
               0, 0, 1,
               0, 0, 1,
               0, 1, 0,
               0, 1, 0), byrow=TRUE, ncol=3)
colnames(mat) = c('A0', 'A1', 'A2')
rownames(mat) = paste0("S", 1:5)
df = as.data.frame(mat)
df$y = c(11, 10.1, 5.2, -0.3, 2)

summary(lm(y~A1+A2, data=df))

## Linear Regression 2 paths - Perfect Linear Relationship with Covariate
mat = matrix(c(1, 0,
               1, 0,
               1, 0,
               0, 1,
               0, 1), byrow=TRUE, ncol=2)
colnames(mat) = c('A0', 'A1')
rownames(mat) = paste0("S", 1:5)
df = as.data.frame(mat)
df$y = c(11, 10.1, 5.2, -0.3, 2)
df$c = c(0.2, 10, 2, 5, 11)

model = lm(y~A1+c, data=df)
summary(model)
lh.o = linearHypothesis(model, c("A1=0"))
lh.o[['Pr(>F)']][2]

## Linear Regression pseudo inversion 3 paths - Moderate with Covariate
mat = matrix(c(1, 0, 0,
               0, 0, 1,
               0, 0, 1,
               0, 1, 0,
               0, 1, 0), byrow=TRUE, ncol=3)
colnames(mat) = c('A0', 'A1', 'A2')
rownames(mat) = paste0("S", 1:5)
df = as.data.frame(mat)
df$y = c(11, 10.1, 5.2, -0.3, 2)
df$c = c(0.2, 10, 2, 5, 11)

model = lm(y~A1+A2+c, data=df)
summary(model)
lh.o = linearHypothesis(model, c("A1=0", "A2=0"))
lh.o
lh.o[['Pr(>F)']][2]

glm(y~A1+A2+c, data=df)

model.n = lm(y~c, data=df)
sse.f = sum(model$residuals^2)
sse.n = sum(model.n$residuals^2)
Fstat = (sse.n - sse.f) / (model.n$df - model$df)
Fstat = Fstat / (sse.f / model$df)
pf(Fstat, model.n$df - model$df, model$df, lower.tail=FALSE)

## logistic regression sandbox
df = data.frame(x=runif(100), y=sample(c(T, F), 100, TRUE))
summary(glm(y~x, family="binomial", data=df))
df$y = rnorm(100, df$x, .1) > .5

df = data.frame(y=sample(c(T, F), 10, TRUE))
df$x = as.numeric(df$y)
N = 1
df$x[1:N] = 1-df$x[1:N]
summary(glm(y~x, family="binomial", data=df))

##
#### Fisher Chi2
##

## Binary phenotype filters minor allele frequency filters
ctab = matrix(c(2, 1,
                6, 0), ncol=2, byrow=TRUE)
fisher.test(ctab)
chisq.test(ctab, simulate.p.value=FALSE, correct=FALSE)

## four alleles and some very low counts
ctab = matrix(c(2, 2, 2, 0,
                1, 0, 2, 1), ncol=4, byrow=TRUE)
chisq.test(ctab, simulate.p.value=FALSE, correct=FALSE)


stoat.group_paths = '33:42,167:158'
ctab = unlist(strsplit(stoat.group_paths, ','))
ctab = lapply(ctab, function(x) as.numeric(unlist(strsplit(x, ':'))))
ctab = do.call(rbind, ctab)

chi2.o = chisq.test(ctab, simulate.p.value=FALSE, correct=FALSE)
chi2.o$p.value

fisher.test(ctab)

ctab = matrix(c(10, 15, 5,
                20, 10, 10), ncol=3)
chisq.test(ctab, simulate.p.value=FALSE, correct=FALSE)
