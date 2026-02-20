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

// #define DEBUG_STATS_TEST

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

// compute the sigmoid function, here corresponding to the "predicted" probability associated to a set of observations and the model (Beta x X)
Eigen::VectorXd LogisticRegression::sigmoid(const Eigen::VectorXd& t) {
    Eigen::VectorXd res = Eigen::VectorXd::Constant(t.size(), 0);
    for (size_t idx = 0; idx < t.size(); idx++) {
        res(idx) = 1 / (1 + std::exp(-t(idx)));
    }
    return res;
}
    
// logistic regression using the Maximum Likelihood Estimate with Newton-Raphson method
test_result_t LogisticRegression::logistic_regression(const BinaryPhenotypeTable& pheno, GenoTable& geno, const CovariateTable& covariates, const double maf, const size_t min_individuals) {
    // prepare an output objet and init to NA
    test_result_t tres;
    tres.pv = "NA";
    // JEAN will the genotype table include samples with no alleles? If yes, we'll need to filter them out (assuming no for now)
    // link the phenotype
    geno.link_to_binary_phenotype(pheno);
    geno.link_to_covariates(covariates);
    // remove non-variable predictors, e.g. alleles absent in all samples
    geno.remove_constant_predictors();    
    geno.remove_noncovered_samples();    
    // should we test this snarl?
    if (!geno.passes_filters(maf, min_individuals)){
        // JEAN what should we do here? returning NA for now
        stoat::LOG_DEBUG("Filtered: didn't pass the filters");
        return tres;
    }
    // add the allele path info to include in the output
    tres.allele_paths = geno.allele_paths_as_str();

    // add covariate with the number of alleles (if necessary) to correct for the parent snarl effect (or normalize)
    geno.add_total_allele_count_covariable();

    // before performing the regression, try to reduce potential colinearity
    geno.remove_duplicated_predictors();    
    geno.remove_one_allele();
    
    // prepare the matrices for the regression
    Eigen::MatrixXd X = geno.make_matrixXd_features();
    Eigen::VectorXd Y = geno.make_vectorxd_phenotype();
#ifdef DEBUG_STATS_TEST
    std::cerr << "X:\n" << X << "\n";
    std::cerr << "Y:\n" << Y << "\n";
#endif
    size_t num_samples = Y.size();
    size_t num_features = X.cols();
    size_t num_predictors = geno.get_n_active_alleles();
    size_t num_covariates = geno.get_n_active_columns() - num_predictors;

    // we'll look for the beta coefficients that maximize the likelihood using the Newton-Raphson method
    // it's an iterative approach so we start by initializing the coefficients (to 0)
    Eigen::VectorXd beta = Eigen::VectorXd::Constant(X.cols(), 0);
#ifdef DEBUG_STATS_TEST
    std::cerr << "beta:\n" << beta << "\n";
#endif

    // should we use Firth's penalized regression? Not at first
    bool penalize = false;

    // prepare matrices that we use multiple times below
    Eigen::MatrixXd Xt = X.transpose();
    // the score, i.e. the gradient or first derivative of the log-likelihood
    // (that's what we are trying to find the root of)
    Eigen::VectorXd score;
    double max_score;
    // how to update the betas, i.e. the delta or step to add to the beta vector
    Eigen::VectorXd delta;
    double max_delta;
    // the predicted probs from the current beta coefficients
    Eigen::VectorXd Ypred;
    // the W matrix with S(BXi)(1-S(BXi)) on the diag
    Eigen::MatrixXd W;
    // the hessian and its inverse
    Eigen::MatrixXd XtWX;
    Eigen::MatrixXd XtWXi;
        
    // we'll iterate to find the MLE with the following parameters
    bool converged = false;
    // the maximum step currently used. Init with default from the class
    double cur_max_step = max_step;
    size_t iter = 0;
    while(iter < max_iterations) {
        // compute the predicted probs from the current beta coefficients
        Ypred = sigmoid(X * beta);
        // W is a diag matrix with S(BXi)(1-S(BXi)) on the diag
        // JEAN might be quite big and used for the diag. Could do better here (sparse matrix? manually multiplying X's columns?)
        W = Eigen::MatrixXd::Constant(Y.size(), Y.size(), 0);
        for (size_t ii = 0; ii < Y.size(); ii++) {
            W(ii, ii) = Ypred(ii) * (1 - Ypred(ii));
        }

        XtWX = Xt * W * X;
        XtWXi = inverse(XtWX);
        // if we want to penalize compute the penalty for the current betas
        if (penalize) {
            Eigen::MatrixXd sqrtW = Eigen::MatrixXd::Constant(Y.size(), Y.size(), 0);
            for (size_t ii = 0; ii < Y.size(); ii++) {
                sqrtW(ii, ii) = std::sqrt(Ypred(ii) * (1 - Ypred(ii)));
            }
            Eigen::MatrixXd hat = (sqrtW * X) * XtWXi * (Xt * sqrtW);
            // the penalization uses the diagonal of this hat matrix as penalty
            // update the Ypred used for the score below with that penalty
            for (size_t ii = 0; ii < Y.size(); ii++) {
                Ypred(ii) = Ypred(ii) - (hat(ii, ii) * (0.5 - Ypred(ii)));
            }
        }

        // compute the score
        score = Xt * (Y - Ypred);

        // update the beta coefficients using the Newton-Raphton formula (XtWX is the hessian)
        // first make sure the update is not too big (usually a bad sign)
        delta = XtWXi * score;
        max_delta = delta.cwiseAbs().maxCoeff();
        if (max_delta > cur_max_step) {
            // reduce the delta so that its largest component is cur_max_step
            delta = cur_max_step * delta / max_delta;
        }
        // ready to update the betas
        beta = beta + delta;

#ifdef DEBUG_STATS_TEST
        std::cerr << "score:\n" << score << "\n\n";
        std::cerr << "beta:\n" << beta << "\n\n";
#endif

        // switch to penalizing if betas get too large (and we were not already penalizing)
        // JEAN any beta larger than 10 is suspicious and unexpected, but maybe there is a better way to decide if we should switch to Firth's regression
        if (!penalize && (beta.maxCoeff() > 10 || beta.maxCoeff() < -10)) {
            penalize = true;
            // start again the iteration process
            iter = 0;
            beta = Eigen::VectorXd::Constant(X.cols(), 0);
#ifdef DEBUG_STATS_TEST
            std::cerr << "switching to penalized Firth regression\n";
            std::cerr << "beta:\n" << beta << "\n\n";
#endif
            continue;
        } else {
            iter++;
        }

        max_score = score.cwiseAbs().maxCoeff();
        // stop if it looks like we've converged (small score and small deltas)
        if (max_score < conv_tol && max_delta < conv_tol) {
            converged = true;
            break;
        }
        // if we're at the last iteration and haven't converged try again with smaller
        // steps from scratch
        if (iter == max_iterations && cur_max_step > 1) {
            iter = 0;
            beta = Eigen::VectorXd::Constant(X.cols(), 0);
            cur_max_step--;
        }
    }

    // if it didn't converge, check if the last iteration looked close to a maximum
    // likelihood (small score and small latest beta update)
    if (!converged) {
        if (max_score < 0.1 && max_delta < 0.1) {
            // not that bad, continue with a warning?
            stoat::LOG_WARN("Logistic regression didn't converge to a great fit (max score coeff " + std::to_string(max_score) + ", max delta " + std::to_string(max_delta) + ") but not too bad. Continuing.");
        } else {
            // too far from a good fit, return NAs.
            return tres;
        }
    }

    // compute log-likelihood for current betas
    Eigen::VectorXd Xbeta = X * beta;
    double loglik = 0;
    for (size_t ii = 0; ii < Y.size(); ii++) {
        loglik += Y(ii) * Xbeta(ii) - std::log(1 + std::exp(Xbeta(ii)));
    }
    if (penalize) {
        // if we're using the Firth penalized regression, the log-likelihood needs to
        // include the penalty, here half the determinant of the Fisher information matrix (Xt W X)
        Ypred = sigmoid(Xbeta);
        W = Eigen::MatrixXd::Constant(Y.size(), Y.size(), 0);
        for (size_t ii = 0; ii < Y.size(); ii++) {
            W(ii, ii) = Ypred(ii) * (1 - Ypred(ii));
        }
        XtWX = X.transpose() * W * X;
        loglik += log(XtWX.determinant()) / 2;
    }

#ifdef DEBUG_STATS_TEST
    std::cerr << "log-likelihood full model:" << loglik << "\n\n";
#endif
    
    // to compute a pvalue, we'll also fit a reduced model without the variables of interest
    // same process as above except we directly use either the standard or penalized
    // approach to match whatever was used for the full model. We also reuse the same cur_max_step
    // keep the intercept and covariates only
    Eigen::MatrixXd X0 = Eigen::MatrixXd::Constant(X.rows(), 1 + num_covariates, 1);
    if (X0.cols() > 1) {
        X0.block(0, 1, X.rows(), num_covariates) = X.block(0, 1 + num_predictors, X.rows(), num_covariates);
    }
    beta = Eigen::VectorXd::Constant(X0.cols(), 0);

#ifdef DEBUG_STATS_TEST
    std::cerr << "X0:\n" << X0 << "\n";
#endif
    
    iter = 0;
    while(iter < max_iterations) {
        Ypred = sigmoid(X0 * beta);
        W = Eigen::MatrixXd::Constant(Y.size(), Y.size(), 0);
        for (size_t ii = 0; ii < Y.size(); ii++) {
            W(ii, ii) = Ypred(ii) * (1 - Ypred(ii));
        }
        XtWX = X0.transpose() * W * X0;
        XtWXi = inverse(XtWX);
        if (penalize) {
            Eigen::MatrixXd sqrtW = Eigen::MatrixXd::Constant(Y.size(), Y.size(), 0);
            for (size_t ii = 0; ii < Y.size(); ii++) {
                sqrtW(ii, ii) = std::sqrt(Ypred(ii) * (1 - Ypred(ii)));
            }
            Eigen::MatrixXd hat = (sqrtW * X0) * XtWXi * (X0.transpose() * sqrtW);
            // the penalization uses the diagonal of this hat matrix as penalty
            // update the Ypred used for the score below with that penalty
            for (size_t ii = 0; ii < Y.size(); ii++) {
                Ypred(ii) = Ypred(ii) - (hat(ii, ii) * (0.5 - Ypred(ii)));
            }
        }
        score = X0.transpose() * (Y - Ypred);
        // first make sure the update is not too big (usually a bad sign)
        delta = XtWXi * score;
        max_delta = delta.cwiseAbs().maxCoeff();
        if (max_delta > cur_max_step) {
            // reduce the delta so that its largest component is cur_max_step
            delta = cur_max_step * delta / max_delta;
        }
        // update the betas
        beta = beta + delta;
        iter++;

        // stop if it looks like we've converged
        if (score.cwiseAbs().maxCoeff() < conv_tol && max_delta < conv_tol) {
            break;
        }
    }
    
#ifdef DEBUG_STATS_TEST
    std::cerr << "reduced model beta:\n" << beta << "\n\n";
#endif
    
    // compute the reduced model
    // compute log-likelihood for reduced model
    Xbeta = X0 * beta;
    double loglik0 = 0;
    for (size_t ii = 0; ii < Y.size(); ii++) {
        loglik0 += Y(ii) * Xbeta(ii) - std::log(1 + std::exp(Xbeta(ii)));
    }
    if (penalize) {
        // if we're using the Firth penalized regression, the log-likelihood needs to
        // include the penalty, here half the determinant of the Fisher information matrix (Xt W X)
        Ypred = sigmoid(Xbeta);
        W = Eigen::MatrixXd::Constant(Y.size(), Y.size(), 0);
        for (size_t ii = 0; ii < Y.size(); ii++) {
            W(ii, ii) = Ypred(ii) * (1 - Ypred(ii));
        }
        XtWX = X0.transpose() * W * X0;
        loglik0 += log(XtWX.determinant()) / 2;
    }
#ifdef DEBUG_STATS_TEST
    std::cout << "log-likelihood reduced model:" << loglik0 << "\n\n";
#endif
    
    // compute -2 * log-likelihood ratio of those models
    double loglik_ratio = -2 * (loglik0 - loglik);
    // should follow a chi2 distribution with df_full-df_reduced degrees of freedom
    int df = X.cols() - X0.cols();

#ifdef DEBUG_STATS_TEST
    std::cerr << "Chi2: " << loglik_ratio << " with df=" << df << "\n\n";
#endif

    if (loglik_ratio < 0) {
        stoat::LOG_WARN("Negative log_likelihood ratio, likely due to issues fitting the full/reduced models. Skipping.");
        return tres;
    }
    
    // compute p-value (double first)
    boost::math::chi_squared_distribution<double> dist(df);
    double p_value = 1.0 - boost::math::cdf(dist, loglik_ratio);

    // try again with more resolution if the pvalue is 0 or infinite (underflow?)
    if (p_value == 0.0 || !std::isfinite(p_value)) {
        boost::multiprecision::cpp_dec_float_50 chi2_hp = loglik_ratio;
        boost::multiprecision::cpp_dec_float_50 df_hp   = df;

        boost::math::chi_squared_distribution<boost::multiprecision::cpp_dec_float_50> dist_hp(df_hp);
        boost::multiprecision::cpp_dec_float_50 p_hp = boost::multiprecision::cpp_dec_float_50(1) - boost::math::cdf(dist_hp, chi2_hp);
        tres.pv = stoat::set_precision_float_50(p_hp);
    } else {
        tres.pv = stoat::set_precision(p_value);
    }

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

test_result_t FisherChi2::fisher_chi2(const BinaryPhenotypeTable& pheno, GenoTable& geno, const double maf, const size_t min_individuals) {
    // prepare an output objet and init to NA
    test_result_t tres;
    tres.pv = "NA";
    tres.second_pv = "NA";
    // JEAN will the genotype table include samples with no alleles? If yes, we'll need to filter them out (assuming no for now)
    // link to the phenotype
    geno.link_to_binary_phenotype(pheno);
    // remove non-variable allele, e.g. absent in both groups
    geno.remove_constant_predictors();    
    geno.remove_noncovered_samples();    
    // should we test this snarl?
    if (!geno.passes_filters(maf, min_individuals)){
        // JEAN what should we do here? returning NA for now
        stoat::LOG_DEBUG("Filtered: didn't pass the filters");
        return tres;
    }
    // fill up the contingency table (one vector per group)
    std::vector<size_t> g0;
    std::vector<size_t> g1;
    geno.fill_contingency_table(g0, g1);

    // performs the test
    auto fc_res = fisher_chi2(g0, g1);
    tres.pv = fc_res.first;
    tres.second_pv = fc_res.second;
    tres.group_paths = stoat::format_group_paths(g0, g1);
    
    return (tres);
}

// ------------------------ Linear regression ------------------------

// Compute Moore-Penrose pseudoinverse using SVD
Eigen::MatrixXd pseudo_inverse(const Eigen::MatrixXd& A, double tol) {
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

// Invert a square matrix
// if the matrix looks non-invertible, use the pseudo-inverse (see above)
Eigen::MatrixXd inverse(const Eigen::MatrixXd &A) {
    if (A.determinant() < 1e-6) {
        return pseudo_inverse(A);
    } else {
        return A.inverse();
    }
}

// Performs linear regression and F-test for predictors only
test_result_t LinearRegression::linear_regression(const QuantitativePhenotypeTable& pheno, GenoTable& geno, const CovariateTable& covariates, const double maf, const size_t min_individuals) {
    // prepare an output objet and init to NA
    test_result_t tres;
    tres.pv = "NA";
    // link the phenotype
    geno.link_to_quantitative_phenotype(pheno);
    geno.link_to_covariates(covariates);
    // remove non-variable allele, e.g. absent in both groups
    geno.remove_constant_predictors();    
    geno.remove_noncovered_samples();    
    // should we test this snarl?
    if (!geno.passes_filters(maf, min_individuals)){
        // JEAN what should we do here? returning NA for now
        stoat::LOG_DEBUG("Filtered: didn't pass the filters");
        return tres;
    }
    // add the allele path info to include in the output
    tres.allele_paths = geno.allele_paths_as_str();

    // add covariate with the number of alleles (if necessary) to correct for the parent snarl effect (or normalize)
    geno.add_total_allele_count_covariable();

    // before performing the regression, try to reduce potential colinearity
    geno.remove_duplicated_predictors();
    geno.remove_one_allele();
    
    Eigen::MatrixXd X_full = geno.make_matrixXd_features();
    Eigen::VectorXd y = geno.make_vectorxd_phenotype();
#ifdef DEBUG_STATS_TEST
    std::cerr << "X:\n" << X_full << "\n";
    std::cerr << "Y:\n" << y << "\n";
#endif
    size_t num_samples = y.size();
    size_t num_predictors = geno.get_n_active_alleles();
    size_t num_covariates = geno.get_n_active_columns() - num_predictors;
    
    stoat::LOG_TRACE("Linear regression. " + std::to_string(num_samples) + " samples, " + std::to_string(num_predictors) + " predictors, " + std::to_string(num_covariates) + " covariates.");
    
    // fit the full model
    // Solve beta X = Y using QR decomposition
    // colPivHouseholderQr is faster than fullPivHouseholderQr but less stable. could switch to full if we encounter problem (or try the HouseolderQR which is less stable but even faster)
    Eigen::VectorXd beta = X_full.colPivHouseholderQr().solve(y);
    Eigen::VectorXd y_hat = X_full * beta;

#ifdef DEBUG_STATS_TEST
    std::cerr << "beta:\n" << beta << "\n";
#endif

    // sum of squared errors and total sum of squares total (like using the mean as prediction)
    double SSE_full = 0.0;
    double SST = 0.0;
    double y_mean = y.sum() / y.size();
    for (int i = 0; i < num_samples; ++i) {
        SSE_full += (y(i) - y_hat(i)) * (y(i) - y_hat(i));
        SST += (y[i] - y_mean) * (y[i] - y_mean);
    }

    // reduced model without any predictor of interest (allele counts)
    double SSE_reduced = 0.0;
    if (num_covariates > 0) {
        // keep the intercept and covariates only
        Eigen::MatrixXd X_reduced = Eigen::MatrixXd::Constant(X_full.rows(), 1 + num_covariates, 1);
        X_reduced.block(0, 1, X_full.rows(), num_covariates) = X_full.block(0, 1 + num_predictors, X_full.rows(), num_covariates);
#ifdef DEBUG_STATS_TEST
        std::cerr << "X reduced:\n" << X_reduced << "\n";
#endif

        Eigen::VectorXd beta_r = X_reduced.colPivHouseholderQr().solve(y);
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
    if (df_denominator <= 0) {
        stoat::LOG_WARN("Too few samples (" + std::to_string(num_samples) + ") compared to alleles+covariates (" + std::to_string(num_params_full) + ") in this snarl. Skipping. Note: increasing the minimum number of individuals with -I could help avoiding those issues and get more robust associations in general.");
        return tres;
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
