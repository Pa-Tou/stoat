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

## read the genotypes
gen.df = read.table(paste0('../../', cur.lab, '/snarl_genotypes.tsv'),
                    comment.char='', header=TRUE, as.is=TRUE, skip=7)
gen.df = gen.df %>% mutate(snarl=paste0(gsub('^.', '', X.START_NODE), '_', gsub('^.', '', END_NODE)))
gen.df = gen.df[,c('snarl', grep('samp', colnames(gen.df), value=TRUE))]

## extract genotypes for a particular snarl and convert to a data.frame compatible with regression
snarl.n = '288_291'
reg.df = gen.df %>% filter(snarl==snarl.n) %>% pivot_longer(-snarl) %>%
  filter(value!='.') %>% 
  mutate(sample=gsub('\\..*', '', name), value=paste0('A', value)) %>%
  group_by(sample, value) %>% summarize(ac=n()) %>%
  pivot_wider(names_from=value, values_from=ac, values_fill=0)

## read phenotype and merge
pheno = read.table('../test_data/input_data/quantitative/phenotype.tsv',
                   header=TRUE, as.is=TRUE)
colnames(pheno) = c('sample', 'pheno')
reg.df = merge(reg.df, pheno)

## should we add a column with the total allele count?
table(reg.df$A0 + reg.df$A1)
reg.df$c = reg.df$A0 + reg.df$A1

## fit a linear regression model
summary(lm(pheno ~ A1 + c, data=reg.df))
