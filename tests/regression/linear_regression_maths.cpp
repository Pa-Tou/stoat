#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono> // for benchmarking
#include <sstream>
#include <string>
#include <map>
#include <numeric>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/distributions/normal.hpp>

// Multiply matrix A (m×n) with matrix B (n×p) -> result is m×p
std::vector<std::vector<double>> matmul(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B) {
    int m = A.size();       // Number of rows in A
    int n = A[0].size();    // Number of columns in A (and rows in B)
    int p = B[0].size();    // Number of columns in B
    // Initialize result matrix with zeros, dimensions m×p
    std::vector<std::vector<double>> result(m, std::vector<double>(p, 0.0));
    
    // Perform matrix multiplication: result[i][j] = sum over k of A[i][k] * B[k][j]
    for (int i = 0; i < m; ++i)
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < p; ++j)
                result[i][j] += A[i][k] * B[k][j];
    return result;
}

// Multiply matrix A (m×n) with vector b (length n) -> result is vector of size m
std::vector<double> matvec(const std::vector<std::vector<double>> &A, const std::vector<double> &b) {
    int m = A.size();       // Number of rows in A
    int n = A[0].size();    // Number of columns in A (should match size of vector b)
    std::vector<double> result(m, 0.0);  // Initialize result vector with zeros
    
    // Compute matrix-vector multiplication: result[i] = sum over j of A[i][j] * b[j]
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            result[i] += A[i][j] * b[j];
    return result;
}

// Compute transpose of a matrix A (m×n) -> result is n×m
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> &A) {
    int m = A.size();       // Number of rows in A
    int n = A[0].size();    // Number of columns in A
    std::vector<std::vector<double>> result(n, std::vector<double>(m, 0.0)); // Transposed matrix n×m
    
    // Assign element result[j][i] = A[i][j]
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            result[j][i] = A[i][j];
    return result;
}

// Invert a square matrix A (n×n) using naive Gaussian elimination without pivoting
std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>> &A) {
    int n = A.size();  // Dimension of square matrix
    std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));  // Identity matrix
    std::vector<std::vector<double>> B = A;  // Copy of matrix A for manipulation
    
    // Initialize identity matrix I
    for (int i = 0; i < n; ++i) {
        I[i][i] = 1.0;
    }

    // Perform Gaussian elimination to convert B to identity and I to inverse of A
    for (int i = 0; i < n; ++i) {
        double diag = B[i][i];
        // Normalize the pivot row by dividing by the pivot element
        for (int j = 0; j < n; ++j) {
            B[i][j] /= diag;
            I[i][j] /= diag;
        }
        // Eliminate all other entries in the current column
        for (int k = 0; k < n; ++k) {
            if (k == i) continue;  // Skip pivot row
            double factor = B[k][i];
            for (int j = 0; j < n; ++j) {
                B[k][j] -= factor * B[i][j];  // Make B[k][i] zero
                I[k][j] -= factor * I[i][j];  // Apply same operations to identity matrix
            }
        }
    }
    // After elimination, I holds the inverse of A
    return I;
}

void linear_regression(
    const std::vector<std::vector<double>>& X_raw,
    const std::vector<double>& y,
    const std::vector<std::vector<double>>& covariates) {

    int n = X_raw.size();
    int num_variants = X_raw[0].size();
    int num_covariates = covariates.empty() ? 0 : covariates[0].size();
    int num_features = 1 + num_variants + num_covariates; // intercept + variants + covariates

    // Build design matrix with intercept, X_raw, and covariates
    std::vector<std::vector<double>> X(n, std::vector<double>(num_features, 1.0));
    for (int i = 0; i < n; ++i) {
        int col = 1;
        for (int j = 0; j < num_variants; ++j)
            X[i][col++] = X_raw[i][j];
        for (int j = 0; j < num_covariates; ++j)
            X[i][col++] = covariates[i][j];
    }

    // OLS computations
    auto Xt = transpose(X);
    auto XtX = matmul(Xt, X);
    auto XtX_inv = inverse(XtX);
    auto Xty = matvec(Xt, y);

    std::vector<double> beta(num_features, 0.0);
    for (int i = 0; i < num_features; ++i)
        for (int j = 0; j < num_features; ++j)
            beta[i] += XtX_inv[i][j] * Xty[j];

    std::vector<double> y_hat = matvec(X, beta);
    double sse = 0.0;
    for (int i = 0; i < n; ++i)
        sse += (y[i] - y_hat[i]) * (y[i] - y_hat[i]);

    // Compute mean of y
    double y_mean = std::accumulate(y.begin(), y.end(), 0.0) / y.size();

    // Compute SST (total sum of squares)
    double sst = 0.0;
    for (int i = 0; i < n; ++i)
        sst += (y[i] - y_mean) * (y[i] - y_mean);

    // Compute R²
    double r2 = (sst == 0.0) ? 1.0 : 1.0 - sse / sst;

    double df_resid = (n - num_features > 0) ? n - num_features : 1;
    double sigma2 = sse / df_resid;

    // --- F-test computation (no p-value, no covariate/intercept) ---
    double ssr = sst - sse;  // Regression sum of squares
    double msr = ssr / num_variants;
    double mse = sse / df_resid;
    double f_stat = msr / mse;
    
    boost::math::fisher_f dist(num_variants, df_resid);
    double p_value = 1 - boost::math::cdf(dist, f_stat);

    std::cout << "p_value : " << p_value << std::endl;
}

int main() {

    std::vector<std::vector<double>> df = {
        {1, 0},
        {1, 0},
        {1, 0},
        {1, 0},
        {1, 0},
        {1, 0},
        {1, 0},
        {0, 1},
        {0, 0},
    };

    std::vector<double> quantitative_phenotype = {4.5, 7.0, 9.2, 10.9, 13.0, 14.0, 11.0, 15.0, 16.0};
    std::vector<std::vector<double>> covariates = {};
    linear_regression(df, quantitative_phenotype, covariates);
    return 0;
}

// g++ -std=c++17 -lboost_math_c99 -o simple_linear linear_regression_maths.cpp
// ./simple_linear

// beta[0] = 16, SE = 8.16561, t = 1.95944, p = 0.300417
// beta[1] = -6.05714, SE = 8.7294, t = -0.693878, p = 0.613823
// beta[2] = -1, SE = 11.5479, t = -0.0865957, p = 0.945009