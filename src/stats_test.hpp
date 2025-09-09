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

#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/distributions/normal.hpp>

#include "arg_parser.hpp"
#include "matrix.hpp"
#include "utils.hpp"

using namespace std;
using boost::multiprecision::cpp_dec_float_50;

namespace stoat {

// ------------------------ Regression class ------------------------

class FisherKhi2 {
    public:
        FisherKhi2(size_t degrees_of_freedom = 1);
        ~FisherKhi2() = default;

        // Function to perform the Chi-square test on row size > 2 
        std::string chi2_2xN(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

        // Function to perform the Chi-square test on row size == 2 
        std::string chi2_2x2(const size_t& m11, const size_t& m12,
            const size_t& m21, const size_t& m22);

        // Function to perform Fisher's exact test
        // not const& because we change the value
        std::string fastFishersExactTest(size_t m11, size_t m12,
            size_t m21, size_t m22);
        
        std::pair<std::string, std::string> fisher_khi2(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

    private:
        // Constants with maximum usable precision for 'double'
        static constexpr double kExactTestEpsilon2 = 9.094947017729282e-13;
        static constexpr double kExactTestBias = 1.0339757656912846e-25;

        // Chi-squared distribution
        const boost::math::chi_squared chi_squared_dist;
        const boost::math::chi_squared_distribution<cpp_dec_float_50> cpp_dec_float_50_dist;
};

class LinearRegression {
    public:
        LinearRegression() = default;
        ~LinearRegression() = default;

        std::tuple<std::string, std::string, std::string, std::string> linear_regression(
            const std::vector<std::vector<double>>& X_raw,
            const std::vector<double>& y,
            const std::vector<std::vector<double>>& covariates);
        
        std::vector<std::vector<double>> inverse(
            const std::vector<std::vector<double>> &A, 
            double tol = 1e-10);

        std::vector<std::vector<double>> pseudoInverse(
            const std::vector<std::vector<double>>& A, 
            double tol = 1e-10);
        
        std::vector<std::vector<double>> fromEigenMatrix(
            const Eigen::MatrixXd& mat);
        
        Eigen::MatrixXd toEigenMatrix(
            const std::vector<std::vector<double>>& mat);
        
        std::vector<std::vector<double>> transpose(
            const std::vector<std::vector<double>> &A);

        std::vector<double> matvec(
            const std::vector<std::vector<double>> &A, 
            const std::vector<double> &b);

        std::vector<std::vector<double>> matmul(
            const std::vector<std::vector<double>> &A, 
            const std::vector<std::vector<double>> &B);
};

class LogisticRegression {
    public:
        LogisticRegression() = default;
        ~LogisticRegression() = default;

        double calculate_log_likelihood(const Eigen::VectorXd& y, const Eigen::VectorXd& p);

        // Standard normal cumulative distribution function
        double normal_cdf(double z);

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

class LMM {
    public:
        LMM() = default;
        ~LMM() = default;

        // template <typename T> 
        // lmm(const std::vector<std::vector<double>>& df,
        //     const std::vector<T>& phenotype_table,
        //     const stoat_vcf::KinshipMatrix& kinship,
        //     const std::vector<std::vector<double>>& covariates);
};

} // namespace stoat

#endif 
