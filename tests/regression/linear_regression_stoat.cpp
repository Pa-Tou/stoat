// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp14)]]

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT

#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <tuple>
#include <iomanip>
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <chrono>

#include <boost/math/distributions/fisher_f.hpp>

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

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

// Convert std::vector<std::vector<double>> to Eigen::MatrixXd
Eigen::MatrixXd toEigenMatrix(const std::vector<std::vector<double>>& mat) {
    int rows = mat.size();
    int cols = mat[0].size();
    Eigen::MatrixXd result(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result(i, j) = mat[i][j];
    return result;
}

// Convert Eigen::MatrixXd back to std::vector<std::vector<double>>
std::vector<std::vector<double>> fromEigenMatrix(const Eigen::MatrixXd& mat) {
    int rows = mat.rows();
    int cols = mat.cols();
    std::vector<std::vector<double>> result(rows, std::vector<double>(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result[i][j] = mat(i, j);
    return result;
}

// Compute Moore-Penrose pseudoinverse using SVD
std::vector<std::vector<double>> pseudoInverse(const std::vector<std::vector<double>>& A, double tol=1e-10) {
    Eigen::MatrixXd mat = toEigenMatrix(A);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);

    const auto& U = svd.matrixU();
    const auto& V = svd.matrixV();
    const auto& S = svd.singularValues();

    Eigen::MatrixXd S_pinv(mat.cols(), mat.rows());
    S_pinv.setZero();

    for (int i = 0; i < S.size(); ++i) {
        if (S(i) > tol)
            S_pinv(i, i) = 1.0 / S(i);
    }

    Eigen::MatrixXd A_pinv = V * S_pinv * U.transpose();
    return fromEigenMatrix(A_pinv);
}

// Invert a square matrix (naive Gaussian elimination, no pivoting)
std::vector<std::vector<double>> inverse(
    const std::vector<std::vector<double>> &A, double tol=1e-10) {

    int n = A.size();
    std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> B = A;

    for (int i = 0; i < n; ++i)
        I[i][i] = 1.0;

    for (int i = 0; i < n; ++i) {
        double diag = B[i][i];

        // Check for near-zero pivot
        if (std::abs(diag) < tol) {
            // Matrix is likely singular or rank-deficient
            return pseudoInverse(A);
        }

        // Normalize the pivot row
        for (int j = 0; j < n; ++j) {
            B[i][j] /= diag;
            I[i][j] /= diag;
        }

        // Eliminate other rows
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

// [[Rcpp::export]]
List cpp_linear_regression_stoat(NumericMatrix Xr, NumericVector yr) {

    // Convert R data to std::vector<std::vector<double>>
    std::vector<std::vector<double>> X_raw(Xr.nrow(), std::vector<double>(Xr.ncol()));
    for (int i = 0; i < Xr.nrow(); ++i)
        for (int j = 0; j < Xr.ncol(); ++j)
            X_raw[i][j] = Xr(i,j);

    std::vector<double> y(yr.begin(), yr.end());

    int num_samples = X_raw.size();
    int num_variants = X_raw[0].size();
    int num_features = 1 + num_variants; // intercept + variants + covariates

    // Build design matrix with intercept, X_raw, and covariates
    std::vector<std::vector<double>> X(num_samples, std::vector<double>(num_features, 1.0));
    for (int i = 0; i < num_samples; ++i) {
        int col = 1;
        for (int j = 0; j < num_variants; ++j)
            X[i][col++] = X_raw[i][j];
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
    for (int i = 0; i < num_samples; ++i)
        sse += (y[i] - y_hat[i]) * (y[i] - y_hat[i]);

    // Compute mean of y
    double y_mean = std::accumulate(y.begin(), y.end(), 0.0) / y.size();

    // Compute SST (total sum of squares)
    double sst = 0.0;
    for (int i = 0; i < num_samples; ++i)
        sst += (y[i] - y_mean) * (y[i] - y_mean);

    // Compute R²
    double r2 = (sst == 0.0) ? 1.0 : 1.0 - sse / sst; // R squared
    double df_resid = (num_samples - num_features - 1 > 0) ? num_samples - num_features - 1 : 1; // Residual Degrees of Freedom
    double sigma2 = sse / df_resid; // Residual Variance

    // --- F-test computation (no p-value, no covariate/intercept) ---
    double ssr = sst - sse;  // Regression sum of squares
    double msr = ssr / num_variants; // Mean Square Regression
    double mse = sse / df_resid; // Mean Squared Error
    double f_stat = msr / mse;  // F-statistic

    boost::math::fisher_f dist(num_variants, df_resid); // F-distribution
    double p_value = 1 - boost::math::cdf(dist, f_stat); // p-value for F-test

    // std::cout << "P-value: " << p_value << std::endl;
    // std::cout << "R²: " << r2 << std::endl;
    // std::cout << "Residual Degrees of Freedom: " << df_resid << std::endl;
    // std::cout << "Mean Squared Error (MSE): " << mse << std::endl;

    return List::create(
        _["r2"] = r2,
        _["p_value"] = p_value
    );

}

// Error in `data.frame()`:
// ! argument manquant, sans valeur associée par défaut
// Backtrace:
//     x
//  1. \-base::data.frame(...)

// Quitting from linear_simulation.rmd:116-147 [setup]
// Exécution arrêtée
