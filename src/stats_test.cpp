#include "stats_test.hpp"

using boost::multiprecision::cpp_dec_float_50;
using boost::math::chi_squared_distribution;

// Fisher's Exact Test for 2x2 contingency table
#ifndef DBL_MAX
#  define DBL_MAX 1.7976931348623157e308
#endif

#ifdef __cplusplus
#  define K_CAST(type, val) (const_cast<type>(val))
#  define R_CAST(type, val) (reinterpret_cast<type>(val))
#  define S_CAST(type, val) (static_cast<type>(val))
#endif

namespace stoat {

// ------------------------ Logistic regression ------------------------

// Standard normal cumulative distribution function
double LogisticRegression::normal_cdf(double z) {
    static const boost::math::normal_distribution<> standard_normal(0.0, 1.0);
    return boost::math::cdf(standard_normal, z);
}

// Sigmoid function
inline double LogisticRegression::sigmoid(double x) {
    return 1.0 / (1.0 + std::exp(-x));
}

// Clamp helper
inline double LogisticRegression::clamp(double x, double lo, double hi) {
    return std::max(lo, std::min(hi, x));
}

// Log-likelihood
double LogisticRegression::calculate_log_likelihood(const Eigen::VectorXd& y, const Eigen::VectorXd& p) {
    double epsilon = 1e-8;
    double ll = 0.0;
    for (int i = 0; i < y.size(); ++i) {
        double pi = clamp(p(i), epsilon, 1.0 - epsilon);
        ll += y(i) * std::log(pi) + (1 - y(i)) * std::log(1 - pi);
    }
    return ll;
}

// GLM Implementation with Iteratively Reweighted Least Squares (IRLS)
std::tuple<std::string, std::string, std::string> LogisticRegression::logistic_regression(
    const std::vector<std::vector<double>>& df,
    const std::vector<bool>& phenotype,
    const std::vector<std::vector<double>>& covariates) {

    size_t num_samples = df.size();
    size_t num_variants = df[0].size();
    size_t num_covariates = 0;
    size_t num_features = num_variants + 1; // +1 for intercept
    
    if (!covariates.empty()) {
        size_t num_covariates = covariates[0].size();
        size_t num_features =  num_variants + num_covariates + 1; // +1 for intercept
    }

    Eigen::MatrixXd X(num_samples, num_features);
    X.col(0) = Eigen::VectorXd::Ones(num_samples);  // Intercept column
    Eigen::VectorXd y(num_samples);
    
    for (size_t i = 0; i < num_samples; ++i) {
        size_t col = 1;

        // Copy variant data
        for (size_t j = 0; j < num_variants; ++j) {
            X(i, col++) = df[i][j];
        }

        for (size_t j = 0; j < num_covariates; ++j) {
            X(i, col++) = covariates[i][j];
        }

        // Binary phenotype
        y(i) = phenotype[i] ? 1.0 : 0.0;
    }

    Eigen::VectorXd beta = Eigen::VectorXd::Zero(num_features);
    Eigen::VectorXd beta_old = beta;
    Eigen::VectorXd p(num_samples);
    Eigen::VectorXd weights(num_samples);

    bool converged = false;
    for (int iter = 0; iter < max_iterations; ++iter) {
        Eigen::VectorXd z = X * beta;
        for (int i = 0; i < num_samples; ++i) {
            p(i) = sigmoid(z(i));
            weights(i) = clamp(p(i) * (1.0 - p(i)), epsilon, 1.0);
        }

        Eigen::MatrixXd X_weighted = X;
        for (int i = 0; i < num_samples; ++i)
            X_weighted.row(i) *= std::sqrt(weights(i));

        Eigen::MatrixXd hessian = X_weighted.transpose() * X_weighted;
        hessian += l2_penalty * Eigen::MatrixXd::Identity(num_features, num_features);

        Eigen::VectorXd gradient = X.transpose() * (y - p) - l2_penalty * beta;

        Eigen::LDLT<Eigen::MatrixXd> ldlt(hessian);
        if (ldlt.info() != Eigen::Success) return std::make_tuple("NA", "NA", "NA");

        Eigen::VectorXd delta = ldlt.solve(gradient);
        beta += delta;

        if ((beta - beta_old).norm() < tolerance) {
            converged = true;
            break;
        }
        beta_old = beta;
    }

    if (!converged) return std::make_tuple("NA", "NA", "NA");

    // Final weights
    Eigen::VectorXd z_final = X * beta;
    for (int i = 0; i < num_samples; ++i) {
        p(i) = sigmoid(z_final(i));
        weights(i) = clamp(p(i) * (1.0 - p(i)), epsilon, 1.0);
    }

    // Covariance matrix
    Eigen::MatrixXd X_weighted = X;
    for (int i = 0; i < num_samples; ++i)
        X_weighted.row(i) *= std::sqrt(weights(i));

    Eigen::MatrixXd hessian = X_weighted.transpose() * X_weighted;
    hessian += l2_penalty * Eigen::MatrixXd::Identity(num_features, num_features);
    Eigen::MatrixXd cov = hessian.inverse();
    Eigen::VectorXd se = cov.diagonal().array().sqrt();

    // --- Wald Test (Normal approximation)
    std::vector<double> p_values;
    p_values.reserve(num_features - 1 - num_covariates);
    for (size_t i = 1; i < num_features - num_covariates; ++i) { // skip intercept
        double z_score = beta(i) / se(i);
        p_values.push_back(2.0 * (1.0 - normal_cdf(std::abs(z_score)))); // Two-sided
    }

    // --- McFadden's R²
    // double ll_full = calculate_log_likelihood(y, p);
    // double p_null_val = clamp(y.mean(), epsilon, 1.0 - epsilon);
    // Eigen::VectorXd p_null = Eigen::VectorXd::Constant(num_samples, p_null_val);
    // double ll_null = calculate_log_likelihood(y, p_null);
    // double r2 = 0.0;
    // if (ll_null != 0.0) {
    //     r2 = 1.0 - (ll_full / ll_null);
    //     r2 = std::min(std::max(r2, 0.0), 1.0);  // clamp between 0 and 1
    // }

    double p_value_adjusted = p_values[0];
    double beta_adjusted = beta(1);
    double se_adjusted = se(1);

    if (p_values.size() > 1) { // case > 3 column/path
        auto [p_values_adjusted, min_index] = stoat::adjusted_hochberg(p_values);
        beta_adjusted = beta[min_index+1];
        se_adjusted = se[min_index+1];
    }

    // set precision : 4 digit
    // std::string r2_str = stoat::set_precision(r2);
    std::string p_value_str = stoat::set_precision(p_value_adjusted);
    std::string beta_str = stoat::set_precision(beta_adjusted);
    std::string se_str = stoat::set_precision(se_adjusted);

    return std::make_tuple(p_value_str, beta_str, se_str);
}

// ------------------------ Chi2 test ------------------------
FisherKhi2::FisherKhi2(size_t degrees_of_freedom) : chi_squared_dist(degrees_of_freedom), cpp_dec_float_50_dist(degrees_of_freedom) {}

std::string FisherKhi2::chi2_2x2(const size_t& a, const size_t& b, const size_t& c, const size_t& d) {

    int64_t row1 = a + b;
    int64_t row2 = c + d;
    int64_t col1 = a + c;
    int64_t col2 = b + d;
    int64_t total = row1 + row2;

    if (row1 == 0 || row2 == 0 || col1 == 0 || col2 == 0) return "NA";

    double expected_a = (double)(row1) * (col1) / total;
    double expected_b = (double)(row1) * (col2) / total;
    double expected_c = (double)(col1) * (row2) / total;
    double expected_d = (double)(col2) * (row2) / total;

    if (expected_a == 0 || expected_b == 0 || expected_c == 0 || expected_d == 0)
        return stoat::set_precision(std::numeric_limits<double>::max());

    double chi2_stat = 0;
    chi2_stat += std::pow((double)a - expected_a, 2) / expected_a;
    chi2_stat += std::pow((double)b - expected_b, 2) / expected_b;
    chi2_stat += std::pow((double)c - expected_c, 2) / expected_c;
    chi2_stat += std::pow((double)d - expected_d, 2) / expected_d;

    if (chi2_stat > 85.0) {
        cpp_dec_float_50 chi2_stat_float_50 = chi2_stat;
        cpp_dec_float_50 pval = 1.0 - boost::math::cdf(cpp_dec_float_50_dist, chi2_stat_float_50);
        return stoat::set_precision_float_50(pval.convert_to<double>());
    }
    return stoat::set_precision(1.0 - boost::math::cdf(chi_squared_dist, chi2_stat));
}

// Check if the observed matrix is valid (no zero rows/columns)
std::string FisherKhi2::chi2_2xN(const std::vector<size_t>& g0, const std::vector<size_t>& g1) {

    size_t cols = g0.size();
    std::vector<size_t> col_totals(cols);
    size_t total = 0;
    size_t row_total_0 = 0;
    size_t row_total_1 = 0;

    for (size_t i = 0; i < cols; ++i) {
        col_totals[i] = g0[i] + g1[i];
        total += col_totals[i];
        row_total_0 += g0[i];
        row_total_1 += g1[i];
    }

    if (total == 0)
        return "NA";
    if (row_total_0 == 0 || row_total_1 == 0)
        return "NA";
    if (std::any_of(col_totals.begin(), col_totals.end(), [](int x){ return x == 0; }))
        return "NA";

    // Compute chi-squared
    double chi2 = 0.0;
    for (size_t i = 0; i < cols; ++i) {
        double expected_0 = static_cast<double>(row_total_0) * col_totals[i] / total;
        double expected_1 = static_cast<double>(row_total_1) * col_totals[i] / total;

        chi2 += (g0[i] - expected_0) * (g0[i] - expected_0) / expected_0;
        chi2 += (g1[i] - expected_1) * (g1[i] - expected_1) / expected_1;
    }

    size_t df = cols - 1;
    if (chi2 > 85.0) { // avoiding case 0.000+00 precision
        cpp_dec_float_50 chi2_stat_float_50 = chi2;
        boost::math::chi_squared_distribution<cpp_dec_float_50> cpp_dec_float_50_dist_2xN(df);
        cpp_dec_float_50 p_value = 1.0 - boost::math::cdf(cpp_dec_float_50_dist_2xN, chi2_stat_float_50);
        return stoat::set_precision_float_50(p_value);
    }

    boost::math::chi_squared dist_2xN(df);
    double pvalue = 1.0 - boost::math::cdf(dist_2xN, chi2);
    return stoat::set_precision(pvalue);
}

// ------------------------ Fisher exact test ------------------------

// Fisher's exact test for a 2x2 contingency table
// m11, m12, m21, m22 are the counts in the table
// Returns the p-value as a std::string with 4 decimal places
std::string FisherKhi2::fastFishersExactTest(size_t m11, size_t m12,
                                 size_t m21, size_t m22) {
    
    // Check for any full-zero row or column
    if ((m11 | m12) == 0 || (m21 | m22) == 0 || (m11 | m21) == 0 || (m12 | m22) == 0) {
        return "NA";
    }
    
    double tprob = (1 - kExactTestEpsilon2) * kExactTestBias;
    double cur_prob = tprob;
    double cprob = 0;
    size_t uii;
    double cur11, cur12, cur21, cur22;
    double preaddp;

    // Ensure we are left of the distribution center, m11 <= m22, and m12 <= m21.
    if (m12 > m21) {
        uii = m12;
        m12 = m21;
        m21 = uii;
    }

    if (m11 > m22) {
        uii = m11;
        m11 = m22;
        m22 = uii;
    }
    
    if ((S_CAST(size_t, m11) * m22) > (S_CAST(size_t, m12) * m21)) {
        uii = m11;
        m11 = m12;
        m12 = uii;
        uii = m21;
        m21 = m22;
        m22 = uii;
    }

    cur11 = m11;
    cur12 = m12;
    cur21 = m21;
    cur22 = m22;

    while (cur12 > 0.5) {
        cur11 += 1;
        cur22 += 1;
        cur_prob *= (cur12 * cur21) / (cur11 * cur22);
        cur12 -= 1;
        cur21 -= 1;
        if (cur_prob > DBL_MAX) {
        return "0";
        }
        if (cur_prob < kExactTestBias) {
        tprob += cur_prob;
        break;
        }
        cprob += cur_prob;
    }

    if (cprob == 0) {
        return "1";
    }

    while (cur12 > 0.5) {
        cur11 += 1;
        cur22 += 1;
        cur_prob *= (cur12 * cur21) / (cur11 * cur22);
        cur12 -= 1;
        cur21 -= 1;
        preaddp = tprob;
        tprob += cur_prob;
        if (tprob <= preaddp) {
        break;
        }
    }

    if (m11) {
        cur11 = m11;
        cur12 = m12;
        cur21 = m21;
        cur22 = m22;
        cur_prob = (1 - kExactTestEpsilon2) * kExactTestBias;
        do {
        cur12 += 1;
        cur21 += 1;
        cur_prob *= (cur11 * cur22) / (cur12 * cur21);
        cur11 -= 1;
        cur22 -= 1;
        preaddp = tprob;
        tprob += cur_prob;
        if (tprob <= preaddp) {
            return stoat::set_precision(preaddp / (cprob + preaddp));
        }
        } while (cur11 > 0.5);
    }

    return stoat::set_precision(tprob / (cprob + tprob));
}

std::pair<std::string, std::string> FisherKhi2::fisher_khi2(const std::vector<size_t>& g0, const std::vector<size_t>& g1) {
    
    std::string chi2_p_value = "NA";
    std::string fastfisher_p_value = "NA";
    
    // Compute  Fisher's exact & Chi-squared test p-value
    if (g0.size() == 2) {
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        chi2_p_value = chi2_2x2(a, b, c, d);
        fastfisher_p_value = fastFishersExactTest(a, b, c, d);
    } else {
        chi2_p_value = chi2_2xN(g0, g1);
    }

    return {chi2_p_value, fastfisher_p_value};
}
// ------------------------ Linear regression ------------------------

Eigen::MatrixXd 
    LinearRegression::pseudoInverse(
    const Eigen::MatrixXd& X, 
    double tol = 1e-6) {

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd singularValues_inv = svd.singularValues();

    for (int i = 0; i < singularValues_inv.size(); ++i) {
        singularValues_inv(i) = (singularValues_inv(i) > tol) ? 1.0 / singularValues_inv(i) : 0.0;
    }

    return svd.matrixV() * singularValues_inv.asDiagonal() * svd.matrixU().transpose();
}

Eigen::MatrixXd 
    LinearRegression::computeXtXinverse(
    const Eigen::MatrixXd& X, 
    double tol = 1e-10) {

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
        LOG_DEBUG("Using pseudo-inverse.");
        return pseudoInverse(XtX);
    }

    // Return actual inverse: XtX⁻¹ = ldlt.solve(I)
    return ldlt.solve(Eigen::MatrixXd::Identity(XtX.rows(), XtX.cols()));
}

// Linear regression function OLS with intercept + covariate
std::tuple<std::string, std::string, std::string, std::string> 
    LinearRegression::linear_regression(
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
        if (std::isnan(t_stats[i]) || std::isinf(t_stats[i])) { // Special case
            LOG_WARN("t-statistics is NaN or Inf, likely due to zero variance in genotype data. Setting p-value to 1.0.");
            p_values.push_back(1.0); // Assign a high p-value for invalid t-statistics
            LOG_DEBUG("Invalid t-statistic encountered");
            continue;
        }
        p_values.push_back(2 * boost::math::cdf(boost::math::complement(t_dist, std::abs(t_stats[i])))); // two-tailed
    }

    double p_value_adjusted = p_values[0];
    double beta_adjusted = beta[1];
    double se_adjusted = se[1];

    if (p_values.size() > 1) {
        auto [p_values_adjusted, min_index] = stoat::adjusted_hochberg(p_values);
        beta_adjusted = beta[min_index+1];
        se_adjusted = se[min_index+1];
    }

    // set precision : 4 digit
    std::string p_value_str = stoat::set_precision(p_value_adjusted);
    std::string beta_str = stoat::set_precision(beta_adjusted);
    std::string se_str = stoat::set_precision(se_adjusted);
    std::string r2_str = stoat::set_precision(r2);

    return std::make_tuple(p_value_str, beta_str, se_str, r2_str);
}

} // namespace stoat