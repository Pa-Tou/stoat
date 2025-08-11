import sys
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

# --- Load input files from command line ---
feature_file = sys.argv[1]
pheno_file = sys.argv[2]

# --- Load features ---
features = pd.read_csv(feature_file, sep='\t')
features = features.rename(columns={features.columns[0]: "IID"})

# --- Load phenotype ---
pheno = pd.read_csv(pheno_file, sep='\t', usecols=["IID", "PHENO"])

# --- Merge on sample ID ---
df = pd.merge(features, pheno, on="IID")

# --- Prepare data for model ---
X = df.drop(columns=["IID", "PHENO"])
# X = sm.add_constant(X)  # Adds intercept term
y = df["PHENO"]

# --- Fit logistic regression model ---
model = sm.Logit(y, X)
result = model.fit(disp=False)

# --- Output coefficients and p-values ---
print(result.summary())  # Shows coef, std err, z, p-value, conf int
print("result.pvalues : ", result.pvalues)

# --- Extract p-values ---
pvals = result.pvalues

# --- Apply Bonferroni correction ---
reject, pvals_corrected, _, _ = multipletests(pvals, method='bonferroni')

# --- Report results ---
print("\nBonferroni-corrected results:")
for coef, pval, pval_corr, rej in zip(pvals.index, pvals, pvals_corrected, reject):
    print(f"{coef}: p = {pval:.4g}, corrected p = {pval_corr:.4g}, significant: {rej}")