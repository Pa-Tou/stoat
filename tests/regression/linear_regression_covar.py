import statsmodels.api as sm
import pandas as pd

# Define your dataset
X_main = [
    [0.5, 0, 0.5],
    [0, 0.5, 0.5],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0.5, 0]
]

covariates = [
    [1.0],
    [2.0],
    [1.0],
    [3.0],
    [2.0]
]

y = [10.5, 13.0, 15.8, 19.7, 21.5]

# Convert to DataFrame
df_main = pd.DataFrame(X_main, columns=['X1', 'X2', 'X3'])
df_cov = pd.DataFrame(covariates, columns=['Cov1'])

# Combine features
X_full = pd.concat([df_main, df_cov], axis=1)

# Add intercept
X_full = sm.add_constant(X_full)

# Fit model
model = sm.OLS(y, X_full).fit()

# Show summary
print(model.summary())
