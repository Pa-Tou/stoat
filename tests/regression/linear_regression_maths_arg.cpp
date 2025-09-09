#include <string>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <sstream>
#include <cmath>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <numeric>

#include <Eigen/Dense>
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
std::vector<std::vector<double>> pseudoInverse(const std::vector<std::vector<double>>& A, double tol = 1e-10) {
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
std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>> &A, double tol = 1e-10) {
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
            // LOG_DEBUG("Using pseudo-inverse.");
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
    boost::math::students_t dist(df_resid);
    std::vector<double> p_values_vector(num_variants, 0.0);
    std::vector<double> beta_vector(num_variants, 0.0);
    std::vector<double> se_vector(num_variants, 0.0);

    for (int i = 1; i < num_variants+1; ++i) {
        double safe_diagonal = XtX_inv[i][i] > 0 ? XtX_inv[i][i] : 0.0;
        double se = std::sqrt(sigma2 * safe_diagonal);
        double t_stat = beta[i] / se;
        double pval;
        if (std::isnan(t_stat) || std::isinf(t_stat)) { // Special case
            pval = 1.0; // Assign a high p-value for invalid t-statistics
            // LOG_DEBUG("Invalid t-statistic encountered");
            continue;
        } else {
            pval = 2 * boost::math::cdf(boost::math::complement(dist, std::fabs(t_stat)));
        }

        // Store results
        beta_vector[i-1] = beta[i];
        se_vector[i-1] = se;
        p_values_vector[i-1] = pval;

        // Print results
        std::cout << "beta[" << i << "] = " << beta[i]
            << ", SE = " << se
            << ", t = " << t_stat
            << ", p = " << pval << '\n';
    }

    std::cout << "R² = " << r2 << std::endl;

    double p_value_adjusted = p_values_vector[0];
    double beta_adjusted = beta_vector[0];
    double se_adjusted = se_vector[0];

    if (p_values_vector.size() > 1) {
        // auto [p_values_adjusted, min_index] = stoat::adjusted_hochberg(p_values);
        // beta_adjusted = beta_vector[min_index];
        // se_adjusted = se_vector[min_index];
    }

    // set precision : 4 digit
    // std::string p_value_str = stoat::set_precision(p_value_adjusted);
    // std::string beta_str = stoat::set_precision(beta_adjusted);
    // std::string se_str = stoat::set_precision(se_adjusted);
    // std::string r2_str = stoat::set_precision(r2);

}

// Function to parse the feature file
void parse_feature_file(
    const std::string& feature_filename,
    std::vector<std::string>& sample_ids,
    std::vector<std::vector<double>>& features) {

    std::ifstream infile(feature_filename);
    if (!infile) {
        throw std::runtime_error("Unable to open feature file");
    }

    std::string line;
    std::getline(infile, line);  // skip header

    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string token;

        std::string sample_id;
        std::getline(ss, sample_id, '\t');
        sample_ids.push_back(sample_id);

        std::vector<double> feature_row;
        while (std::getline(ss, token, '\t')) {
            feature_row.push_back(std::stod(token));
        }

        features.push_back(feature_row);
    }
}

// Function to parse the phenotype file
void parse_phenotype_file(
    const std::string& phenotype_filename,
    const std::vector<std::string>& sample_ids,
    std::vector<double>& phenotype) {

    std::ifstream infile(phenotype_filename);
    if (!infile) {
        throw std::runtime_error("Unable to open phenotype file");
    }

    std::unordered_map<std::string, double> phenotype_map;

    std::string line;
    std::getline(infile, line);  // skip header

    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string fid, iid, pheno_str;
        std::getline(ss, fid, '\t');
        std::getline(ss, iid, '\t');
        std::getline(ss, pheno_str, '\t');

        phenotype_map[iid] = std::stod(pheno_str);
    }

    for (const auto& sample : sample_ids) {
        if (phenotype_map.find(sample) != phenotype_map.end()) {
            phenotype.push_back(phenotype_map[sample]);
        } else {
            throw std::runtime_error("Sample ID not found in phenotype file: " + sample);
        }
    }
}

// Example usage
int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <feature_file> <phenotype_file>\n";
        return 1;
    }

    std::string feature_file = argv[1];
    std::string phenotype_file = argv[2];

    std::vector<std::string> sample_ids;
    std::vector<std::vector<double>> features;
    std::vector<double> phenotype;
    std::vector<std::vector<double>> covar;

    try {
        parse_feature_file(feature_file, sample_ids, features);
        parse_phenotype_file(phenotype_file, sample_ids, phenotype);

        std::cout << "Parsed " << features.size() << " samples with " << features[0].size() << " features.\n";
        std::cout << "Parsed " << phenotype.size() << " phenotype values.\n";

        linear_regression(features, phenotype, covar);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return EXIT_SUCCESS;
}

// g++ -std=c++17 -I/usr/include/eigen3 -lboost_math_c99 -o lr_maths_arg linear_regression_maths_arg.cpp
// ./lr_maths_arg ../../output_droso/regression/5124567_5124564.tsv ../../../lab/droso/data/pangenome_pheno.tsv
