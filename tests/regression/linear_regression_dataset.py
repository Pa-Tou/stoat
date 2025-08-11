import numpy as np
import statsmodels.api as sm

def linear_regression_summary(df, quantitative_phenotype):
    X = np.array(df)
    y = np.array(quantitative_phenotype)

    # Add intercept
    X = sm.add_constant(X)

    try:
        model = sm.OLS(y, X)
        results = model.fit()

        # Return full summary text
        return results.summary()

    except Exception as e:
        return f"Error during regression: {e}"

# Test cases data
test_cases = [
    {
        "name": "Linear Regression 1 - Perfect Linear Relationship",
        "df": [
            [0],
            [1],
            [0]
        ],
        "quantitative_phenotype": [2.0, 4.0, 6.0]
    },
    {
        "name": "Linear Regression 2 - Moderate",
        "df": [
            [0.5, 0],
            [0, 0.5],
            [1, 0],
            [0, 1],
            [0, 0.5]
        ],
        "quantitative_phenotype": [10.5, 13.0, 15.8, 19.7, 21.5]
    },
    {
        "name": "Linear Regression 3 - Weaker Correlation",
        "df": [
                [1, 0],
                [1, 0],
                [1, 0],
                [1, 0],
                [1, 0],
                [1, 0],
                [1, 0],
                [0, 1],
                [0, 0]
        ],
        "quantitative_phenotype": [4.5, 7.0, 9.2, 10.9, 13.0, 14.0, 11.0, 15.0, 16.0]
    }
]

# Run and print summary for each test
for test in test_cases:
    print(f"Test: {test['name']}")
    summary = linear_regression_summary(test["df"], test["quantitative_phenotype"])
    print(summary)
    print("="*80)

# python3 linear_regression_dataset.py

#                             OLS Regression Results                            
# ==============================================================================
# Dep. Variable:                      y   R-squared:                       0.421
# Model:                            OLS   Adj. R-squared:                  0.228
# Method:                 Least Squares   F-statistic:                     2.184
# Date:                Wed, 06 Aug 2025   Prob (F-statistic):              0.194
# Time:                        16:15:46   Log-Likelihood:                -21.782
# No. Observations:                   9   AIC:                             49.56
# Df Residuals:                       6   BIC:                             50.16
# Df Model:                           2                                         
# Covariance Type:            nonrobust                                         
# ==============================================================================
#                  coef    std err          t      P>|t|      [0.025      0.975]
# ------------------------------------------------------------------------------
# const         16.0000      3.334      4.800      0.003       7.843      24.157
# x1            -6.0571      3.564     -1.700      0.140     -14.777       2.663
# x2            -1.0000      4.714     -0.212      0.839     -12.536      10.536
# ==============================================================================
# Omnibus:                        1.016   Durbin-Watson:                   0.443
# Prob(Omnibus):                  0.602   Jarque-Bera (JB):                0.410
# Skew:                          -0.494   Prob(JB):                        0.815
# Kurtosis:                       2.659   Cond. No.                         7.22
# ==============================================================================