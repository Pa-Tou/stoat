## We have some "*system*" tests where we run STOAT on larger files and compare the results.
## However, when we change something major in the code, it might change the pvalues a bit, for example. 
## This script helped double-check that the changes in pvalues were minor, or zoom in to the cases where it changed enough to warrant some investigations.

library(dplyr)
library(ggplot2)
library(tidyr)

## read the expected output
exp.lab = 'output_binary_covar'
exp.df = read.table(paste0('../test_data/expected_output/vcf/', exp.lab, '/stoat.assoc.pvalues.tsv'),
                    comment.char='', header=TRUE, as.is=TRUE)
## read the current output (by manually running the command in the test log)
cur.lab = 'output_binary'
cur.df = read.table(paste0('../../', cur.lab, '/stoat.assoc.pvalues.tsv'),
                    comment.char='', header=TRUE, as.is=TRUE)
head(cur.df)

## merge both
eqtl = TRUE
if(eqtl) {
  df = merge(exp.df, cur.df, by=c('SNARL', 'GENE'), suffixes=c('.e', '.c'))
} else {
  df = merge(exp.df, cur.df, by=c('SNARL'), suffixes=c('.e', '.c'))
}
head(df)


summary(df$P.e-df$P.c)

ggplot(df, aes(x=P.e, P.c)) + geom_point(alpha=.4) + theme_bw()

ggplot(df, aes(x=P.e, P.c)) + geom_point(alpha=.4) + theme_bw() +
  scale_x_log10() + scale_y_log10()

## top variants with most different pvalues
df %>% arrange(desc(abs(P.e-P.c))) %>% head

##
## Investigate the differences
##
library(logistf) ## to test Firth regression as an more robust logistic regression
library(car)
library(lmtest)

## read the genotypes
gen.df = read.table(paste0('../../', cur.lab, '/snarl_genotypes.tsv'),
                    comment.char='', header=TRUE, as.is=TRUE, skip=7)
gen.df = gen.df %>% mutate(snarl=paste0(gsub('^.', '', X.START_NODE), '_', gsub('^.', '', END_NODE)))
gen.df = gen.df[,c('snarl', grep('samp', colnames(gen.df), value=TRUE))]

## extract genotypes for a particular snarl and convert to a data.frame compatible with regression
snarl.n = '650_647'
reg.df = gen.df %>% filter(snarl==snarl.n) %>% pivot_longer(-snarl) %>%
  filter(value!='.') %>% 
  mutate(sample=gsub('\\..*', '', name), value=paste0('A', value)) %>%
  group_by(sample, value) %>% summarize(ac=n()) %>%
  pivot_wider(names_from=value, values_from=ac, values_fill=0)

## read phenotype and merge
pheno = read.table('../test_data/input_data/binary/phenotype.tsv',
                   header=TRUE, as.is=TRUE)
colnames(pheno) = c('sample', 'pheno')

reg.df = merge(reg.df, pheno)

## should we add a column with the total allele count?
table(reg.df$A0 + reg.df$A1)
reg.df$c = reg.df$A0 + reg.df$A1

## fit a linear regression model
summary(lm(pheno ~ A1 + c, data=reg.df))

summary(glm(pheno ~ A1 + c, data=reg.df, family="binomial"))

firth = TRUE
mod.f = logistf(pheno~A1+c, data=reg.df, firth=firth)
mod.f2 = logistf(pheno~A1, data=reg.df, firth=firth)
mod.n = logistf(pheno~1, data=reg.df, firth=firth)
mod.n2 = logistf(pheno~c, data=reg.df, firth=firth)

pchisq(-2*(mod.n2$loglik[1]-mod.f$loglik[1]), 1, lower.tail=FALSE)
pchisq(-2*(mod.n$loglik[1]-mod.f$loglik[1]), 2, lower.tail=FALSE)
pchisq(-2*(mod.n$loglik[1]-mod.f2$loglik[1]), 2, lower.tail=FALSE)

mod.f = glm(pheno ~ A1 + c, data=reg.df, family="binomial")
mod.f2 = glm(pheno ~ A1, data=reg.df, family="binomial")
mod.n = glm(pheno ~ 1, data=reg.df, family="binomial")
mod.n2 = glm(pheno ~ c, data=reg.df, family="binomial")

lrtest(mod.f, mod.n2)
lrtest(mod.f, mod.n)

lrtest(mod.f2, mod.n)


## read the covariates
covar = read.table('../test_data/input_data/binary/covariate.tsv',
                   header=TRUE, as.is=TRUE)
colnames(covar)[1] = 'sample'

reg.df = merge(reg.df, covar)

firth = TRUE
mod.f = logistf(pheno~A1+c+PC1+SEX+PC3, data=reg.df, firth=firth)
mod.f2 = logistf(pheno~A1+PC1+SEX+PC3, data=reg.df, firth=firth)
mod.n = logistf(pheno~PC1+SEX+PC3, data=reg.df, firth=firth)
mod.n2 = logistf(pheno~c+PC1+SEX+PC3, data=reg.df, firth=firth)

pchisq(-2*(mod.n2$loglik[1]-mod.f$loglik[1]), 1, lower.tail=FALSE)
pchisq(-2*(mod.n$loglik[1]-mod.f$loglik[1]), 2, lower.tail=FALSE)
pchisq(-2*(mod.n$loglik[1]-mod.f2$loglik[1]), 2, lower.tail=FALSE)

mod.f = glm(pheno ~ A1 + c+PC1+SEX+PC3, data=reg.df, family="binomial")
mod.f2 = glm(pheno ~ A1+PC1+SEX+PC3, data=reg.df, family="binomial")
mod.n = glm(pheno ~ PC1+SEX+PC3, data=reg.df, family="binomial")
mod.n2 = glm(pheno ~ c+PC1+SEX+PC3, data=reg.df, family="binomial")

lrtest(mod.f, mod.n2)
lrtest(mod.f, mod.n)

lrtest(mod.f2, mod.n)

