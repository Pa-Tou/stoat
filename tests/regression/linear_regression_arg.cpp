#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>  // For t-distribution
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

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Linear regression function OLS with intercept + covariate
void linear_regression(
    const std::vector<std::vector<double>>& df,
    const std::vector<double>& quantitative_phenotype,
    const std::vector<std::vector<double>>& covar) {

    size_t num_samples = df.size();
    size_t num_variants = df[0].size();
    size_t num_covariates = 0;
    size_t num_features = num_variants + 1; // +1 for intercept

    if (!covar.empty()) {
        size_t num_covariates = covar[0].size();
        size_t num_features = num_variants + num_covariates + 1; // +1 for intercept
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
    
    // Coefficients beta
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
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
    Eigen::MatrixXd cov_matrix = (X.transpose() * X).inverse();    
    Eigen::VectorXd se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();

    // change cov_matrix calcul if X.transpose() * X might be ill-conditioned or nearly singular
    if (se.hasNaN()) {
        Eigen::MatrixXd XtX = X.transpose() * X;
        Eigen::MatrixXd cov_matrix = XtX.ldlt().solve(Eigen::MatrixXd::Identity(X.cols(), X.cols()));
        se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();
    }

    // t-statistics
    Eigen::VectorXd t_stats = beta.array() / se.array();
    boost::math::students_t t_dist(df_res);
 
    std::vector<double> p_values;
    for (int i = 0; i < num_features; ++i) { // i = 1 avoid const p-value
        if (std::isnan(t_stats[i]) || std::isinf(t_stats[i])) {
            p_values.push_back(1.0); // Assign a high p-value for invalid t-statistics
            continue;
        }
        p_values.push_back(2 * boost::math::cdf(boost::math::complement(t_dist, std::abs(t_stats[i])))); // two-tailed
        cout << "p_values[" << i << "] : " << p_values[i] << std::endl;
    }

    // Print results
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Coefficients (beta):" << std::endl;
    for (int i = 0; i < num_features; ++i) {
        std::cout << "beta[" << i << "] = " << beta[i] << std::endl;
    }
    std::cout << "Standard Errors (se):" << std::endl;
    for (int i = 0; i < num_features; ++i) {
        std::cout << "se[" << i << "] = " << se[i] << std::endl;
    }
    std::cout << "R²: " << r2 << std::endl;
    std::cout << "Residual Degrees of Freedom: " << df_res << std::endl;
    std::cout << "Mean Squared Error (MSE): " << mse << std::endl;
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

// LINUX
// g++ -std=c++17 -I/usr/include/eigen3 -lboost_math_c99 -o lr_arg linear_regression_arg.cpp

// MACOS
// g++ -std=c++17 -I/usr/local/eigen3 -lboost_math_c99 -o lr_arg linear_regression_arg.cpp

// ./lr_arg ../../output/regression/4_6.tsv ../../data/quantitative/phenotype.tsv
