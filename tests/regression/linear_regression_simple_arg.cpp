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

    double df_resid = (n - p > 0) ? n - p : 1;
    double sigma2 = sse / df_resid;
    
    for (int i = 0; i < p; ++i) {
        double safe_diagonal = XtX_inv[i][i] > 0 ? XtX_inv[i][i] : 0.0;
        double se = std::sqrt(sigma2 * safe_diagonal);
        double t_stat = beta[i] / se;
        boost::math::students_t dist(df_resid);
        double pval = 2 * boost::math::cdf(boost::math::complement(dist, std::fabs(t_stat)));

        std::cout << "beta[" << i << "] = " << beta[i]
            << ", SE = " << se
            << ", t = " << t_stat
            << ", p = " << pval << '\n';
    }
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

// g++ -std=c++11 -O2 -lboost_math_c99 -o lr_simple_arg linear_regression_simple_arg.cpp
// ./lr_simple_arg ../../output/regression/4_6.tsv ../../data/quantitative/phenotype.tsv