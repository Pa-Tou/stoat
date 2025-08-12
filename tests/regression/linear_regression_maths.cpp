#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono> // for benchmarking
#include <Rcpp.h>

#include <boost/math/distributions/students_t.hpp>

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

    // WARNING: If X has multicollinearity (one or more columns are linearly dependent),
    //          XtX becomes singular or nearly singular. In that case, inverse(XtX)
    //          may fail (division by zero) or produce numerically unstable results,
    //          because the naive Gaussian elimination here does not use pivoting.
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

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List cpp_linear_regression_maths(Rcpp::NumericMatrix X_raw, Rcpp::NumericVector y) {
    int n = X_raw.nrow();
    int k = X_raw.ncol();
    int p = 1 + k;

    std::vector<std::vector<double>> X(n, std::vector<double>(p, 1.0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < k; ++j)
            X[i][j+1] = X_raw(i,j);

    auto Xt = transpose(X);
    auto XtX = matmul(Xt, X);

    auto XtX_inv = inverse(XtX);
    std::vector<double> y_vec(y.begin(), y.end());
    auto Xty = matvec(Xt, y_vec);

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

    std::vector<double> p_values;

    for (int i = 0; i < p; ++i) {
        double se = std::sqrt(sigma2 * XtX_inv[i][i]);
        double t_stat = beta[i] / se;
        boost::math::students_t dist(df_resid);
        double pval;
        if (std::isnan(t_stat) || std::isinf(t_stat)) {
            pval = 1.0;
        } else {
            try {
                pval = 2 * boost::math::cdf(boost::math::complement(dist, std::abs(t_stat)));
            } catch (...) {
                pval = 1.0; // fallback if Boost still fails
            }
        }

        p_values.push_back(pval);
    }

    return Rcpp::List::create(
        _["coefficients"] = beta,
        _["p_values"] = p_values
    );
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