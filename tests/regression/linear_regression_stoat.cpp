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
#include <boost/multiprecision/cpp_dec_float.hpp>

using boost::multiprecision::cpp_dec_float_50;

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
List cpp_linear_regression_stoat(NumericMatrix Xr, NumericMatrix Xcovar, NumericVector yr) {

    int n = Xr.nrow();
    int p = Xr.ncol();
    int q = Xcovar.ncol();

    std::vector<double> y(yr.begin(), yr.end());

    // ---- FULL MODEL ----
    int k_full = 1 + p + q;
    std::vector<std::vector<double>> Xfull(n, std::vector<double>(k_full, 1.0));

    // Add X
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < p; ++j)
            Xfull[i][1 + j] = Xr(i, j);

    // Add covariates (if any)
    if (q > 0) {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < q; ++j)
                Xfull[i][1 + p + j] = Xcovar(i, j);
    }

    // ---- Fit FULL model ----
    auto Xt = transpose(Xfull);
    auto XtX = matmul(Xt, Xfull);
    auto XtXi = inverse(XtX);
    auto Xty = matvec(Xt, y);
    auto beta = matvec(XtXi, Xty);
    auto yhat = matvec(Xfull, beta);

    double SSE_full = 0.0;
    for (int i = 0; i < n; ++i)
        SSE_full += (y[i] - yhat[i]) * (y[i] - yhat[i]);

    // ---- Reduced SSE ----
    double SSE_reduced = 0.0;
    if (q > 0) {
        // build reduced model: intercept + covariates
        int k_red = 1 + q;
        std::vector<std::vector<double>> Xred(n, std::vector<double>(k_red, 1.0));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < q; ++j)
                Xred[i][1 + j] = Xcovar(i, j);

        auto Xt_r = transpose(Xred);
        auto XtX_r = matmul(Xt_r, Xred);
        auto XtXi_r = inverse(XtX_r);
        auto Xty_r = matvec(Xt_r, y);
        auto beta_r = matvec(XtXi_r, Xty_r);
        auto yhat_r = matvec(Xred, beta_r);

        for (int i = 0; i < n; ++i)
            SSE_reduced += (y[i] - yhat_r[i]) * (y[i] - yhat_r[i]);
    } else {
        // no covariates: reduced model = intercept only
        double y_mean = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
        for (int i = 0; i < n; ++i)
            SSE_reduced += (y[i] - y_mean) * (y[i] - y_mean);
    }

    // ---- F-statistic ----
    int df_num = p;
    int df_den = n - k_full;
    double num = (SSE_reduced - SSE_full) / df_num;
    double den = SSE_full / df_den;
    double Fstat = num / den;

    cpp_dec_float_50 Fstat_hp = Fstat;
    cpp_dec_float_50 df_num_hp = df_num;
    cpp_dec_float_50 df_den_hp = df_den;

    boost::math::fisher_f_distribution<cpp_dec_float_50> dist_hp(df_num_hp, df_den_hp);
    cpp_dec_float_50 p_value_hp = cpp_dec_float_50(1) - boost::math::cdf(dist_hp, Fstat_hp);

    return List::create(
    _["F_stat"]  = Fstat,
    _["p_value"] = static_cast<double>(p_value_hp),
    _["df_num"]  = df_num,
    _["df_den"]  = df_den
    );
}
