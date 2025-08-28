import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
from statsmodels.stats.outliers_influence import variance_inflation_factor

# Read X (genotype matrix)
X = pd.read_csv("table_test.tsv", sep="\t", header=None)
# remove last column if it's not a feature
# X = X.iloc[:, :-1]

# Read Y (phenotype vector)
Y = pd.read_csv("phenotype.tsv", header=None).squeeze()  # .squeeze() to convert to Series

# Optional: check shapes
if X.shape[0] != Y.shape[0]:
    raise ValueError("X and Y must have the same number of rows")

# Add intercept (required for statsmodels unless already included)
X_with_const = sm.add_constant(X)

# Fit model (optional, just for reference)
model = sm.OLS(Y, X_with_const).fit()
print(model.summary())

# Calculate VIF for each predictor (including constant)
vif_data = pd.DataFrame()
vif_data["feature"] = X_with_const.columns
vif_data["VIF"] = [variance_inflation_factor(X_with_const.values, i)
                   for i in range(X_with_const.shape[1])]

print("\nVIF values:")
print(vif_data)

# Plot VIF
plt.figure(figsize=(10, 5))
plt.bar(vif_data["feature"].astype(str), vif_data["VIF"])
plt.xlabel("Feature")
plt.ylabel("VIF")
plt.title("Variance Inflation Factor (VIF)")
plt.legend()
plt.tight_layout()
plt.show()
