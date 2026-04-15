#ifndef stats_test_HPP
#define stats_test_HPP

#include <string>
#include <vector>
#include <atomic>
#include <Eigen/Dense>
#include <Eigen/Core>

#include "utils.hpp"

namespace stoat {

// used to pass the result of a test (Fisher-Chi2 or regression)
// up to the writer, hence why they're all stings
// maybe not the best way but better than std::pairs and such
// (second_pv can be used for Fisher-Chi2 for the Chi2 pvalue)
struct test_result_t {
    double pv;
    double second_pv;
    std::string group_paths;
    std::string allele_paths;

    std::string to_string(){
        return set_precision(pv) + " " + set_precision(second_pv) + " " + group_paths + " " + allele_paths;
    };
};

// ------------------------ Regression class ------------------------

class FisherChi2 {
public:
    FisherChi2() = default;
    ~FisherChi2() = default;
    
    // Function to perform the Chi-square test on row size > 2 
    // If one of the genotypes has no counts, return "NA"
    // If one of the columns (alleles) has no counts, ignore it
    double chi2_2xN(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

    // Function to perform the Chi-square test on row size == 2 
    double chi2_2x2(const size_t& m11, const size_t& m12,
                         const size_t& m21, const size_t& m22);

    // Function to perform Fisher's exact test
    // not const& because we change the value
    // This is implemented in fast_fishers_exact_test.cpp because it was not us who wrote it
    double fastFishersExactTest(size_t m11, size_t m12,
                                     size_t m21, size_t m22);

    // depending on the number of alleles, uses Chi2 (>2 alleles) or the fast
    // Fisher exact test (2 alleles).
    // returns two pvalues: fastfisher_p_value (which can be NA if >2 alleles), and chi2_p_value
    std::pair<double, double> fisher_chi2(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

private:
    // Constants with maximum usable precision for 'double'
    static constexpr double kExactTestEpsilon2 = 9.094947017729282e-13;
    static constexpr double kExactTestBias = 1.0339757656912846e-25;

    // Atomic: these counters are incremented from parallel OMP threads
    std::atomic<size_t> chi2_zero{0};
    std::atomic<size_t> chi2_inf{0};
};

class LinearRegression {
    public:
        LinearRegression() = default;
        ~LinearRegression() = default;

        double linear_regression(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const size_t num_predictors);

    private:
        // Atomic: these counters are incremented from parallel OMP threads
        std::atomic<size_t> count_number_sse_null{0};
        std::atomic<size_t> count_number_few_sample{0};
        std::atomic<size_t> count_number_f_stat_close_to_0{0};
        std::atomic<size_t> count_number_f_stat_negative{0};
};

class LogisticRegression {
public:
    LogisticRegression() = default;
    ~LogisticRegression() = default;
    
    // compute the sigmoid function, here corresponding to the "predicted" probability associated to a set of observations and the model (Beta x X)
    Eigen::VectorXd sigmoid(const Eigen::VectorXd& t) const;
    
    // Logistic regression using the Newton-Raphson method to find the MLE
    // if the betas get too high, it will switch to Firth penalized regression
    // p-value is computed from the likelihood ratio of the full vs reduced model
    // returns the pvalue as a string, ready to be written to the output
    double logistic_regression(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const size_t num_predictors);
    
private:
    // maximum number of iterations to perform
    const int max_iterations = 100;
    // how big can the step/delta be when updating the betas
    const double max_step = 3;
    // tolerance to decide if the score and delta are small enough to consider the iteration to have converged
    const double conv_tol = 0.001;

    // Atomic: these counters are incremented from parallel OMP threads
    std::atomic<size_t> count_number_not_convergence{0};
    std::atomic<size_t> count_number_negative_log_likelihood{0};

};

// Compute the inverse or Moore-Penrose pseudoinverse of a matrix
Eigen::MatrixXd inverse(const Eigen::MatrixXd& A);
Eigen::MatrixXd pseudo_inverse(const Eigen::MatrixXd& A, double tol = 1e-10);
    
} // namespace stoat

#endif 
