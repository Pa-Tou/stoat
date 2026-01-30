#include "stats_test.hpp"

// using boost::multiprecision::cpp_dec_float_50;
// using boost::math::chi_squared_distribution;

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

// ------------------------ Functions to filter tables ------------------------

// Should we exclude this test? "filter" this table out
bool filter_binary_table(
    std::vector<size_t>& g0, 
    std::vector<size_t>& g1,
    const size_t& individuals_included, 
    const size_t& min_individuals,
    const double& maf_threshold) {

    // not enough individuals/haplotypes OR number of paths < 2
    if (individuals_included < min_individuals || g0.size() < 2) {
        stoat::LOG_DEBUG("Filtered: not enough individuals: " + std::to_string(individuals_included));
        return true; // Empty or invalid input → filter
    }

    // prepare the total allele count
    size_t total_allele_count = 0;
    for (size_t i = 0 ; i < g0.size() ; i++) {
        total_allele_count += g0[i];
        total_allele_count += g1[i];
    }

    // now check if allele frequencies are above the threshold
    int alleles_above_threshold = 0;
    for (size_t i = 0; i < g0.size(); ++i) {
        size_t allele_count = g0[i] + g1[i];
        double allele_freq = static_cast<double>(allele_count) / total_allele_count;
        allele_freq = std::min(allele_freq, 1.0 - allele_freq);

        if (allele_freq > maf_threshold) {
            ++alleles_above_threshold;
        }
    }

    if (alleles_above_threshold < 2) {
        stoat::LOG_DEBUG("Filtered: less than two alleles with frequency above " + std::to_string(maf_threshold));
    }
    // filter if there is not at least two path with frequency > MAF threshold
    return alleles_above_threshold < 2; 
}

// ------------------------ Logistic regression ------------------------

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
test_result_t LogisticRegression::logistic_regression(const BinaryPhenotypeTable& pheno, const GenotypeTable& geno, const CovariateTable& covar, const double maf, const size_t min_individuals) {
    // prepare an output objet and init to NA
    test_result_t tres;
    tres.pv = "NA";
    tres.beta = "NA";
    tres.se = "NA";
    // JEAN will the genotype table include samples with no alleles? If yes, we'll need to filter them out (assuming no for now)
    // combine the phenotype and genotype information
    CombinedTable combined_table(geno);
    combined_table.combine_binary_phenotype(pheno);
    combined_table.combine_covariates(covar);
    // JEAN the functions below could be gathered in a "prepare_for_regression" function, for example, to avoid repetition in the linea regression function
    // remove non-variable predictors, e.g. alleles absent in all samples
    combined_table.remove_constant_predictors();    
    // should we test this snarl?
    if (!combined_table.passes_filters(maf, min_individuals)){
        // JEAN what should we do here? returning NA for now
        stoat::LOG_DEBUG("Filtered: didn't pass the filters");
        return tres;
    }
    // add the allele path info to include in the output
    tres.allele_paths = geno.allele_paths_as_str();

    // add covariate with the number of alleles (if necessary) to correct for the parent snarl effect (or normalize)
    combined_table.add_total_allele_count_covariable();

    // before performing the regression, try to reduce potential colinearity
    combined_table.remove_duplicated_predictors();    
    combined_table.remove_one_allele();
    
    // prepare the matrices for the regression
    Eigen::MatrixXd X = combined_table.make_matrixXd_features();
    Eigen::VectorXd y = combined_table.make_vectorxd_phenotype();
    size_t num_samples = y.size();
    size_t num_features = X.cols();
    size_t num_covariates = combined_table.get_n_covariates();

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

        Eigen::LDLT<Eigen::MatrixXd> ldlt(hessian);
        if (ldlt.info() != Eigen::Success) return tres;

        Eigen::VectorXd gradient = X.transpose() * (y - p) - l2_penalty * beta;
        Eigen::VectorXd delta = ldlt.solve(gradient);
        beta += delta;

        if ((beta - beta_old).norm() < tolerance) {
            converged = true;
            break;
        }
        beta_old = beta;
    }

    if (!converged) return tres;

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
    for (size_t i = 1; i < num_features - num_covariates; ++i) {

        const double z_score = beta(i) / se(i);

        // --- Fast path (double precision) ---
        const boost::math::normal_distribution<double> standard_normal(0.0, 1.0);
        double p_value = 2.0 * (1.0 - boost::math::cdf(standard_normal, std::fabs(z_score)));

        // --- Handle very significant (underflow) cases ---
        if (p_value == 0.0 || !std::isfinite(p_value)) {

            const boost::multiprecision::cpp_dec_float_50 z_hp = std::fabs(z_score);
            const boost::math::normal_distribution<boost::multiprecision::cpp_dec_float_50> standard_normal_hp(0.0, 1.0);
            boost::multiprecision::cpp_dec_float_50 p_hp = boost::multiprecision::cpp_dec_float_50(2) * (boost::multiprecision::cpp_dec_float_50(1) - boost::math::cdf(standard_normal_hp, z_hp));

            // convert to double (or keep string, depending on your output format)
            p_values.push_back(p_hp.convert_to<double>());
        } else {
            p_values.push_back(p_value);
        }
    }
    
    double p_value_adjusted = p_values[0];
    double beta_adjusted = beta(1);
    double se_adjusted = se(1);

    // JEAN this needs to be updated to the "new" F-test, no? (then remove adjusted_hochberg)
    if (p_values.size() > 1) { // case > 3 column/path
        auto [p_values_adjusted, min_index] = stoat::adjusted_hochberg(p_values);
        beta_adjusted = beta[min_index+1];
        se_adjusted = se[min_index+1];
    }

    // set precision : 4 digit
    // std::string r2_str = stoat::set_precision(r2);
    tres.pv = stoat::set_precision(p_value_adjusted);
    tres.beta = stoat::set_precision(beta_adjusted);
    tres.se = stoat::set_precision(se_adjusted);

    return tres;
}

// ------------------------ Chi2 test ------------------------

std::string FisherChi2::chi2_2x2(const size_t& a, const size_t& b, const size_t& c, const size_t& d) {

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

    // Degrees of freedom for 2×2 table
    int df = 1;

    // --- Compute p-value (double first) ---
    boost::math::chi_squared_distribution<double> dist(df);
    double p_value = 1.0 - boost::math::cdf(dist, chi2_stat);

    if (p_value == 0.0 || !std::isfinite(p_value)) {

        // Recompute in high precision if underflow occurs
        boost::multiprecision::cpp_dec_float_50 chi2_hp = chi2_stat;
        boost::multiprecision::cpp_dec_float_50 df_hp   = df;

        boost::math::chi_squared_distribution<boost::multiprecision::cpp_dec_float_50> dist_hp(df_hp);
        boost::multiprecision::cpp_dec_float_50 p_hp = boost::multiprecision::cpp_dec_float_50(1) - boost::math::cdf(dist_hp, chi2_hp);
        return stoat::set_precision_float_50(p_hp);
    }

    return stoat::set_precision(p_value);
}

// Check if the observed matrix is valid (no zero rows/columns)
std::string FisherChi2::chi2_2xN(const std::vector<size_t>& g0, const std::vector<size_t>& g1) {

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

    // Compute chi-squared
    double chi2 = 0.0;
    size_t skipped_cols_count = 0;
    for (size_t i = 0; i < cols; ++i) {

        // Skip any allele with no sample counts
        if (col_totals[i] == 0) { 
            skipped_cols_count++;
            continue;
        }
        double expected_0 = static_cast<double>(row_total_0) * col_totals[i] / total;
        double expected_1 = static_cast<double>(row_total_1) * col_totals[i] / total;

        chi2 += (g0[i] - expected_0) * (g0[i] - expected_0) / expected_0;
        chi2 += (g1[i] - expected_1) * (g1[i] - expected_1) / expected_1;
    }

    size_t df = cols - skipped_cols_count - 1;

    // --- Fast path: double precision ---
    boost::math::chi_squared_distribution<double> dist(df);
    double p_value = 1.0 - boost::math::cdf(dist, chi2);

    // If underflow or invalid, recompute in high precision
    if (p_value == 0.0 || !std::isfinite(p_value)) {

        boost::multiprecision::cpp_dec_float_50 chi2_hp = chi2;
        boost::multiprecision::cpp_dec_float_50 df_hp   = df;
        boost::math::chi_squared_distribution<boost::multiprecision::cpp_dec_float_50> dist_hp(df_hp);

        boost::multiprecision::cpp_dec_float_50 p_hp = boost::multiprecision::cpp_dec_float_50(1) - boost::math::cdf(dist_hp, chi2_hp);
        return stoat::set_precision_float_50(p_hp);
    }

    return stoat::set_precision(p_value);
}

// ------------------------ Fisher exact test ------------------------


std::pair<std::string, std::string> FisherChi2::fisher_chi2(const std::vector<size_t>& g0, const std::vector<size_t>& g1) {
    
    std::string chi2_p_value = "NA";
    std::string fastfisher_p_value = "NA";
    
    // compute  Fisher's exact or Chi-squared test p-value
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

    return {fastfisher_p_value, chi2_p_value};
}

test_result_t FisherChi2::fisher_chi2(const BinaryPhenotypeTable& pheno, const GenotypeTable& geno, const double maf, const size_t min_individuals) {
    // prepare an output objet and init to NA
    test_result_t tres;
    tres.pv = "NA";
    tres.second_pv = "NA";
    // JEAN will the genotype table include samples with no alleles? If yes, we'll need to filter them out (assuming no for now)
    // combine the phenotype and genotype information
    CombinedTable combined_table(geno);
    combined_table.combine_binary_phenotype(pheno);
    // remove non-variable allele, e.g. absent in both groups
    combined_table.remove_constant_predictors();    
    // should we test this snarl?
    if (!combined_table.passes_filters(maf, min_individuals)){
        // JEAN what should we do here? returning NA for now
        stoat::LOG_DEBUG("Filtered: didn't pass the filters");
        return tres;
    }
    // fill up the contingency table (one vector per group)
    int n_alleles = combined_table.get_n_alleles();
    std::vector<size_t> g0(n_alleles, 0);
    std::vector<size_t> g1(n_alleles, 0);
    auto phenotype_vec = combined_table.get_phenotype();
    auto genotype_mat = combined_table.get_predictors();
    assert(genotype_mat.size() == n_alleles); // no covariates allowed here
    // loop through the alleles and samples and update the contingency table
    for (int al_ii=0; al_ii < n_alleles; al_ii++) {
        for (int samp_ii=0; samp_ii < phenotype_vec.size(); samp_ii++) {
            if (genotype_mat[al_ii][samp_ii] > 0){
                // tally the allele counts in each group
                if (phenotype_vec[samp_ii] == 0){
                    g0[al_ii] += genotype_mat[al_ii][samp_ii];
                } else {
                    g1[al_ii] += genotype_mat[al_ii][samp_ii];
                }
            }
        }
    }

    // performs the test
    auto fc_res = fisher_chi2(g0, g1);
    tres.pv = fc_res.first;
    tres.second_pv = fc_res.second;
    tres.group_paths = stoat::format_group_paths(g0, g1);
    
    return (tres);
}

// ------------------------ Linear regression ------------------------

// Compute Moore-Penrose pseudoinverse using SVD
Eigen::MatrixXd LinearRegression::pseudo_inverse(const Eigen::MatrixXd& A, double tol) {
    // SVD
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const auto& U = svd.matrixU();
    const auto& V = svd.matrixV();
    const auto& S = svd.singularValues();

    Eigen::MatrixXd S_pinv(A.cols(), A.rows());
    S_pinv.setZero();

    for (int i = 0; i < S.size(); ++i) {
        if (S(i) > tol)
            S_pinv(i, i) = 1.0 / S(i);
    }

    return V * S_pinv * U.transpose();
}

// Invert a square matrix (naive Gaussian elimination, no pivoting)
// more info that might have inspired this code (?) https://math.uww.edu/~mcfarlat/inverse.htm and https://cppscripts.com/cpp-matrix-inversion
Eigen::MatrixXd LinearRegression::inverse(const Eigen::MatrixXd &A, double tol) {
    size_t nrows = A.rows();

    // the "augmented" matrix
    Eigen::MatrixXd B = A;

    // prepare the identity matrix
    Eigen::MatrixXd I = Eigen::MatrixXd::Constant(nrows, nrows, 0);
    for (int i = 0; i < nrows; ++i) {
        I(i, i) = 1.0;
    }

    // normalize and reduce each pivot (?)
    for (int i = 0; i < nrows; ++i) {
        double diag = B(i, i);

        // Check for near-zero pivot
        if (std::abs(diag) < tol) {
            // Matrix is likely singular or rank-deficient
            LOG_DEBUG("Using pseudo-inverse.");
            return pseudo_inverse(A);
        }

        // Normalize the pivot row
        for (int j = 0; j < nrows; ++j) {
            B(i, j) /= diag;
            I(i, j) /= diag;
        }

        // Eliminate other rows
        for (int k = 0; k < nrows; ++k) {
            if (k == i) continue;
            double factor = B(k, i);
            for (int j = 0; j < nrows; ++j) {
                B(k, j) -= factor * B(i, j);
                I(k, j) -= factor * I(i, j);
            }
        }
    }

    return I;
}

// Performs linear regression and F-test for predictors only
test_result_t LinearRegression::linear_regression(const QuantitativePhenotypeTable& pheno, const GenotypeTable& geno, const CovariateTable& covariates, const double maf, const size_t min_individuals) {
    // prepare an output objet and init to NA
    test_result_t tres;
    tres.pv = "NA";
    tres.r2 = "NA";
    // JEAN will the GenotypeTable include samples with no alleles? If yes, we'll need to filter them out (assuming no for now)
    // combine the phenotype and genotype information
    CombinedTable combined_table(geno);
    combined_table.combine_quantitative_phenotype(pheno);
    combined_table.combine_covariates(covariates);
    // remove non-variable allele, e.g. absent in both groups
    combined_table.remove_constant_predictors();    
    // should we test this snarl?
    if (!combined_table.passes_filters(maf, min_individuals)){
        // JEAN what should we do here? returning NA for now
        stoat::LOG_DEBUG("Filtered: didn't pass the filters");
        return tres;
    }
    // add the allele path info to include in the output
    tres.allele_paths = geno.allele_paths_as_str();

    // add covariate with the number of alleles (if necessary) to correct for the parent snarl effect (or normalize)
    combined_table.add_total_allele_count_covariable();

    // before performing the regression, try to reduce potential colinearity
    combined_table.remove_duplicated_predictors();    
    combined_table.remove_one_allele();
    
    Eigen::MatrixXd X_full = combined_table.make_matrixXd_features();
    Eigen::VectorXd y = combined_table.make_vectorxd_phenotype();
    size_t num_samples = y.size();
    size_t num_predictors = combined_table.get_n_alleles();
    size_t num_covariates = combined_table.get_n_covariates();
    
    stoat::LOG_TRACE("Linear regression. " + std::to_string(num_samples) + " samples, " + std::to_string(num_predictors) + " predictors, " + std::to_string(num_covariates) + " covariates.");
    
    // fit the full model
    Eigen::MatrixXd Xt = X_full.transpose();
    Eigen::MatrixXd XtXi = inverse(Xt * X_full);
    Eigen::VectorXd Xty = Xt * y;
    Eigen::VectorXd beta = XtXi * Xty;
    Eigen::VectorXd y_hat = X_full * beta;

    // sum of squared errors and total sum of squares total (like using the mean as prediction)
    double SSE_full = 0.0;
    double SST = 0.0;
    double y_mean = y.sum() / y.size();
    for (int i = 0; i < num_samples; ++i) {
        SSE_full += (y(i) - y_hat(i)) * (y(i) - y_hat(i));
        SST += (y[i] - y_mean) * (y[i] - y_mean);
    }
    // save the R2
    tres.r2 = stoat::set_precision(1.0 - SSE_full / SST);

    // reduced model without any predictor of interest (allele counts)
    double SSE_reduced = 0.0;
    if (num_covariates > 0) {
        // keep the intercept and covariates only
        Eigen::MatrixXd X_reduced = Eigen::MatrixXd::Constant(X_full.rows(), 1 + num_covariates, 1);
        X_reduced.block(0, 1, X_full.rows(), num_covariates) = X_full.block(0, 1 + num_predictors, X_full.rows(), num_covariates);

        Eigen::MatrixXd Xt_r = X_reduced.transpose();
        Eigen::MatrixXd XtXi_r = inverse(Xt_r * X_reduced);
        Eigen::VectorXd Xty_r = Xt_r * y;
        Eigen::VectorXd beta_r = XtXi_r * Xty_r;
        Eigen::VectorXd y_hat_r = X_reduced * beta_r;

        for (int i = 0; i < num_samples; ++i) {
            SSE_reduced += (y(i) - y_hat_r(i)) * (y(i) - y_hat_r(i));
        }
    } else {
        // no covariates: reduced model = intercept only
        SSE_reduced = SST;
    }

    // time to compute the F-statistic
    double F_stat = 0;
    // F = [(SSE_reduced - SSE_full) / df_numerator] / [SSE_full / df_denominator]
    // Numerator df = number of tested predictors
    int df_numerator = num_predictors;
    // Denominator df = residual df in full model
    int num_params_full = 1 + num_predictors + num_covariates; // intercept + predictors + covariates
    int df_denominator = num_samples - num_params_full;
    
    // JEAN problem can arise if we have less samples than variables
    // maybe we should skip those tests when they happen? For now, warning the user and recommending increasing -I
    if (df_denominator < 0) {
        stoat::LOG_WARN("Less samples (" + std::to_string(num_samples) + ") than alleles+covariates (" + std::to_string(num_params_full) + ") in this snarl. Increasing the minimum number of individuals with -I to get more robust associations and avoid issues.");
    } else {
        double numerator = (SSE_reduced - SSE_full) / df_numerator;
        double denominator = SSE_full / df_denominator;
        F_stat = numerator / denominator;

        if (F_stat < 0 && F_stat > -0.00001){
            stoat::LOG_WARN("F statistic is negative but very close to 0 (" + std::to_string(F_stat) + "). Assuming it's 0.");
            F_stat = 0;
        }
        if (F_stat < -0.00001){
            stoat::LOG_WARN("F statistic is negative: " + std::to_string(F_stat) + ". This is concerning. Recommendation: increase the minimum number of individuals with -I to get more robust associations and avoid issues.");
            F_stat = 0;
        }
    }
    
    // compute a P-value
    boost::math::fisher_f_distribution<double> dist(df_numerator, df_denominator);
    double p_value = 1.0 - boost::math::cdf(dist, F_stat);

    // recompute P-value with higher precision if underflow or invalid result
    if (p_value == 0.0 || !std::isfinite(p_value)) {
        const boost::multiprecision::cpp_dec_float_50 F_hp = F_stat;
        const boost::multiprecision::cpp_dec_float_50 df_n_hp = df_numerator;
        const boost::multiprecision::cpp_dec_float_50 df_d_hp = df_denominator;
        boost::math::fisher_f_distribution<boost::multiprecision::cpp_dec_float_50> dist_hp(df_n_hp, df_d_hp);
        boost::multiprecision::cpp_dec_float_50 p_hp = boost::multiprecision::cpp_dec_float_50(1) - boost::math::cdf(dist_hp, F_hp);
        // convert to string
        tres.pv = stoat::set_precision_float_50(p_hp);
    } else {
        // convert to string
        tres.pv = stoat::set_precision(p_value);
    }
    
    return tres;
}

} // namespace stoat
