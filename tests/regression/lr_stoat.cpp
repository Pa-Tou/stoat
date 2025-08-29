#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <boost/math/distributions/students_t.hpp>
#include <chrono>
#include <Eigen/SVD>

Eigen::MatrixXd pseudoInverse(const Eigen::MatrixXd& X, double tol = 1e-6) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd singularValues_inv = svd.singularValues();

    for (int i = 0; i < singularValues_inv.size(); ++i) {
        singularValues_inv(i) = (singularValues_inv(i) > tol) ? 1.0 / singularValues_inv(i) : 0.0;
    }

    return svd.matrixV() * singularValues_inv.asDiagonal() * svd.matrixU().transpose();
}

Eigen::MatrixXd computeXtXinverse(const Eigen::MatrixXd& X, double tol = 1e-10) {
    Eigen::MatrixXd XtX = X.transpose() * X;
    Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX);

    Eigen::VectorXd D = ldlt.vectorD();
    bool rank_deficient = false;
    for (int i = 0; i < D.size(); ++i) {
        if (std::abs(D(i)) < tol) {
            rank_deficient = true;
            break;
        }
    }

    if (rank_deficient) {
        std::cerr << "Warning: Matrix XᵀX may be rank-deficient. Multicollinearity detected.\n";
        return pseudoInverse(XtX);
    }

    // Return actual inverse: XtX⁻¹ = ldlt.solve(I)
    return ldlt.solve(Eigen::MatrixXd::Identity(XtX.rows(), XtX.cols()));
}

void linear_regression(
    const std::vector<std::vector<double>>& df,
    const std::vector<double>& quantitative_phenotype,
    const std::vector<std::vector<double>>& covar) {

    size_t num_samples = df.size();
    size_t num_variants = df[0].size();
    size_t num_covariates = 0;
    size_t num_features = num_variants + 1; // +1 for intercept

    if (!covar.empty()) {
        num_covariates = covar[0].size();
        num_features = num_variants + num_covariates + 1; // +1 for intercept
    }

    Eigen::MatrixXd X(num_samples, num_features);
    X.col(0) = Eigen::VectorXd::Ones(num_samples);  // Intercept column
    Eigen::VectorXd y(num_samples);
    
    for (size_t i = 0; i < num_samples; ++i) {
        y(i) = quantitative_phenotype[i];
        size_t col = 1;
        for (size_t j = 0; j < num_variants; ++j) {
            X(i, col++) = df[i][j];
        }
        for (size_t j = 0; j < num_covariates; ++j) {
            X(i, col++) = covar[i][j];
        }
    }

    Eigen::MatrixXd XtXinverse = computeXtXinverse(X);

    // Compute beta using the pseudo-inverse
    Eigen::VectorXd beta = XtXinverse * (X.transpose() * y);
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;

    // R²
    double rss = residuals.squaredNorm();
    double tss = (y.array() - y.mean()).matrix().squaredNorm();
    double r2 = 1 - (rss / tss);

    int df_res = (num_samples - X.cols() + 1); // residual degrees of freedom
    df_res = std::max(df_res, 1); // Ensure df_res is at least 1 to avoid division by zero
    double mse = rss / df_res;

    // Standard errors
    Eigen::VectorXd se = (XtXinverse.diagonal() * mse).array().sqrt().matrix();

    // t-statistics
    Eigen::VectorXd t_stats = beta.array() / se.array();
    boost::math::students_t t_dist(df_res);

    std::vector<double> p_values;
    for (int i = 1; i < num_variants+1; ++i) { // i = 1 avoid const p-value
        p_values.push_back(2 * boost::math::cdf(boost::math::complement(t_dist, std::abs(t_stats[i])))); // two-tailed
        std::cout << "p_values[" << i-1 << "] : " << p_values[i-1] << std::endl;
    }

    // Print results
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Coefficients (beta):" << std::endl;
    for (int i = 1; i < num_variants+1; ++i) {
        std::cout << "beta[" << i << "] = " << beta[i] << std::endl;
    }

    std::cout << "Standard Errors (se):" << std::endl;
    for (int i = 1; i < num_variants+1; ++i) {
        std::cout << "se[" << i << "] = " << se[i] << std::endl;
    }

    std::cout << "R²: " << r2 << std::endl;
    std::cout << "Residual Degrees of Freedom: " << df_res << std::endl;
    std::cout << "Mean Squared Error (MSE): " << mse << std::endl;
}

// === MAIN with Example Data ===
int main() {

    std::vector<std::vector<double>> df = {
        {0, 0.5, 0.5},
        {0.5, 0, 0},
        {0.25, 0.5, 0},
        {0.5, 0, 0},
        {0.333333, 0.333333, 0}
    };

    std::vector<double> quantitative_phenotype = {39.8333, 30.24, 49, 72, 54.3673};

    std::vector<std::vector<double>> covariates = {};

    linear_regression(df, quantitative_phenotype, covariates);
    return EXIT_SUCCESS;
}

// LINUX
// g++ -std=c++17 -I/usr/include/eigen3 -lboost_math_c99 -o lr_stoat lr_stoat.cpp

// MACOS
// g++ -std=c++17  -I/usr/local/include/eigen3 -lboost_math_c99 -o lr_stoat lr_stoat.cpp

// ./lr_stoat
