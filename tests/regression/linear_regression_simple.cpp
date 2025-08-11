#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono> // for benchmarking

#include <boost/math/distributions/students_t.hpp>

// Multiply matrix A (m×n) with matrix B (n×p) -> result is m×p
std::vector<std::vector<double>> matmul(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B) {
    int m = A.size(), n = A[0].size(), p = B[0].size();
    std::vector<std::vector<double>> result(m, std::vector<double>(p, 0.0));
    for (int i = 0; i < m; ++i)
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < p; ++j)
                result[i][j] += A[i][k] * B[k][j];
    return result;
}

// Multiply matrix A (m×n) with vector b (n) -> result is vector of size m
std::vector<double> matvec(const std::vector<std::vector<double>> &A, const std::vector<double> &b) {
    int m = A.size(), n = A[0].size();
    std::vector<double> result(m, 0.0);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            result[i] += A[i][j] * b[j];
    return result;
}

// Transpose of a matrix
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> &A) {
    int m = A.size(), n = A[0].size();
    std::vector<std::vector<double>> result(n, std::vector<double>(m, 0.0));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            result[j][i] = A[i][j];
    return result;
}

// Invert a square matrix (naive Gaussian elimination, no pivoting)
std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>> &A) {
    int n = A.size();
    std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> B = A;

    for (int i = 0; i < n; ++i)
        I[i][i] = 1.0;

    for (int i = 0; i < n; ++i) {
        double diag = B[i][i];
        for (int j = 0; j < n; ++j) {
            B[i][j] /= diag;
            I[i][j] /= diag;
        }
        for (int k = 0; k < n; ++k) {
            if (k == i) continue;
            double factor = B[k][i];
            for (int j = 0; j < n; ++j) {
                B[k][j] -= factor * B[i][j];
                I[k][j] -= factor * I[i][j];
            }
        }
    }
    return I;
}

void linear_regression(
    const std::vector<std::vector<double>>& X_raw,
    const std::vector<double>& y,
    const std::vector<std::vector<double>>& covariates) {

    int n = X_raw.size();
    int k = X_raw[0].size();
    int c = covariates.empty() ? 0 : covariates[0].size();
    int p = 1 + k + c; // intercept + variants + covariates

    // Build design matrix with intercept, X_raw, and covariates
    std::vector<std::vector<double>> X(n, std::vector<double>(p, 1.0));
    for (int i = 0; i < n; ++i) {
        int col = 1;
        for (int j = 0; j < k; ++j)
            X[i][col++] = X_raw[i][j];
        for (int j = 0; j < c; ++j)
            X[i][col++] = covariates[i][j];
    }

    // OLS computations
    auto Xt = transpose(X);
    auto XtX = matmul(Xt, X);
    auto XtX_inv = inverse(XtX);
    auto Xty = matvec(Xt, y);

    std::vector<double> beta(p, 0.0);
    for (int i = 0; i < p; ++i)
        for (int j = 0; j < p; ++j)
            beta[i] += XtX_inv[i][j] * Xty[j];

    std::vector<double> y_hat = matvec(X, beta);
    double sse = 0.0;
    for (int i = 0; i < n; ++i)
        sse += (y[i] - y_hat[i]) * (y[i] - y_hat[i]);

    double df_resid = (n - p) <= 0 ? n - p : 1; // avoid Degrees of freedom <= 0
    double sigma2 = sse / df_resid;

    for (int i = 0; i < p; ++i) {
        double se = std::sqrt(sigma2 * XtX_inv[i][i]);
        double t_stat = beta[i] / se;
        boost::math::students_t dist(df_resid);
        double pval = 2 * boost::math::cdf(boost::math::complement(dist, std::fabs(t_stat)));

        std::cout << "beta[" << i << "] = " << beta[i]
            << ", SE = " << se
            << ", t = " << t_stat
            << ", p = " << pval << '\n';
    }
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

// g++ -std=c++17 -lboost_math_c99 -o simple_linear linear_regression_simple.cpp
// ./simple_linear

// beta[0] = 16, SE = 8.16561, t = 1.95944, p = 0.300417
// beta[1] = -6.05714, SE = 8.7294, t = -0.693878, p = 0.613823
// beta[2] = -1, SE = 11.5479, t = -0.0865957, p = 0.945009