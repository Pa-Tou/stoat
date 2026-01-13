#ifndef stats_test_HPP
#define stats_test_HPP

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
#include <Eigen/Dense>
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <cassert>

#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/distributions/normal.hpp>

#include "arg_parser.hpp"
#include "matrix.hpp"
#include "utils.hpp"
#include "feature_tables.hpp"

namespace stoat {

// ------------------------ Regression class ------------------------

class FisherChi2 {
public:
    FisherChi2() = default;
    ~FisherChi2() = default;
    
    // Function to perform the Chi-square test on row size > 2 
    // If one of the genotypes has no counts, return "NA"
    // If one of the columns (alleles) has no counts, ignore it
    std::string chi2_2xN(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

    // Function to perform the Chi-square test on row size == 2 
    std::string chi2_2x2(const size_t& m11, const size_t& m12,
                         const size_t& m21, const size_t& m22);

    // Function to perform Fisher's exact test
    // not const& because we change the value
    // This is implemented in fast_fishers_exact_test.cpp because it was not us who wrote it
    std::string fastFishersExactTest(size_t m11, size_t m12,
                                     size_t m21, size_t m22);

    // depending on the number of alleles, uses Chi2 (>2 alleles) or the fast
    // Fisher exact test (2 alleles).
    // returns two pvalues: chi2_p_value, fastfisher_p_value (which can be NA if >2 alleles)
    // JEAN rename function or flip the returned pvalues to match the order in the function name
    std::pair<std::string, std::string> fisher_chi2(const std::vector<size_t>& g0, const std::vector<size_t>& g1);
    // JEAN don't really like that we are passing filtering thresholds as arguments. Should we make a well-defined struct with input parameters?
    std::pair<std::string, std::string> fisher_chi2(const BinaryPhenotypeTable& pheno, const GenotypeTable& geno, const size_t maf, const size_t min_individuals);

private:
    // Constants with maximum usable precision for 'double'
    static constexpr double kExactTestEpsilon2 = 9.094947017729282e-13;
    static constexpr double kExactTestBias = 1.0339757656912846e-25;
};

class LinearRegression {
    public:
        LinearRegression() = default;
        ~LinearRegression() = default;

        std::tuple<std::string, std::string> linear_regression(
            const std::vector<std::vector<double>>& X_raw,
            const std::vector<double>& y,
            const std::vector<std::vector<double>>& covariates);
        
        std::vector<std::vector<double>> inverse(
            const std::vector<std::vector<double>> &A, 
            double tol = 1e-10);

        std::vector<std::vector<double>> pseudoInverse(
            const std::vector<std::vector<double>>& A, 
            double tol = 1e-10);

        std::vector<std::vector<double>> transpose(
            const std::vector<std::vector<double>> &A);

        std::vector<double> mult_mat_vec(
            const std::vector<std::vector<double>> &A, 
            const std::vector<double> &b);

        std::vector<std::vector<double>> mult_mat_mat(
            const std::vector<std::vector<double>> &A, 
            const std::vector<std::vector<double>> &B);
};

class LogisticRegression {
    public:
        LogisticRegression() = default;
        ~LogisticRegression() = default;

        double calculate_log_likelihood(const Eigen::VectorXd& y, const Eigen::VectorXd& p);

        // Sigmoid function
        inline double sigmoid(double x);

        // Clamp helper
        inline double clamp(double x, double lo, double hi);

        // GLM Implementation with Iteratively Reweighted Least Squares (IRLS)
        std::tuple<std::string, std::string, std::string> logistic_regression(
            const std::vector<std::vector<double>>& variant_data,
            const std::vector<bool>& phenotype,
            const std::vector<std::vector<double>>& covariates);

    private:
        const int max_iterations = 100;
        const double tolerance = 1e-6;
        const double l2_penalty = 1e-4;
        const double epsilon = 1e-8;
};

void combine_identical_columns_quantitative_table(
    std::vector<std::vector<double>>& df);

// Modify g0 and g1 to remove any column (allele/path) where both g0 and g1 are empty
// Return a vector of bools with true if the column was kept and false removed
std::vector<bool> remove_empty_columns_binary_table(
    std::vector<size_t>& g0, 
    std::vector<size_t>& g1);

void remove_empty_columns_quantitative_table(
    std::vector<std::vector<double>>& df);

void remove_last_columns_quantitative_table(
    std::vector<std::vector<double>>& df);

bool check_last_columns_quantitative_table(
    const std::vector<std::vector<double>>& df);

/// Return true if snarl must be filtered
bool filter_quantitative_table(
    const std::vector<std::vector<double>>& df,
    const size_t& min_individuals,
    const double& maf);

/// Return true if the snarl should be filtered out, false if it should be kept
bool filter_binary_table(
    std::vector<size_t>& g0, 
    std::vector<size_t>& g1,
    const size_t& individuals_included,
    const size_t& min_individuals,
    const double& maf);

} // namespace stoat

#endif 
