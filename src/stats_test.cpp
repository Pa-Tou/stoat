#include "stats_test.hpp"
#include "utils.hpp"
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/chi_squared.hpp>

using boost::math::chi_squared_distribution;

// #define DEBUG_STATS_TEST

namespace stoat {

// ------------------------ Logistic regression ------------------------

// compute the sigmoid function, here corresponding to the "predicted" probability associated to a set of observations and the model (Beta x X)
Eigen::VectorXd LogisticRegression::sigmoid(const Eigen::VectorXd& t) const {
    Eigen::VectorXd res = Eigen::VectorXd::Constant(t.size(), 0);
    for (size_t idx = 0; idx < t.size(); idx++) {
        res(idx) = 1 / (1 + std::exp(-t(idx)));
    }
    return res;
}
    
// logistic regression using the Maximum Likelihood Estimate with Newton-Raphson method
double LogisticRegression::logistic_regression(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const size_t num_predictors) const {

#ifdef DEBUG_STATS_TEST
    std::cerr << "X:\n" << X << "\n";
    std::cerr << "Y:\n" << Y << "\n";
#endif
    size_t num_samples = Y.size();
    size_t num_features = X.cols();
    size_t num_covariates = num_features - num_predictors - 1;

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
            stoat::LOG_WARN("Logistic regression didn't converge to a great fit (max score coeff " + 
                std::to_string(max_score) + ", max delta " + std::to_string(max_delta) +
                ") but not too bad. Continuing.", "count_number_not_convergence");
#ifdef DEBUG_STATS_TEST
    std::stringstream ss_X;
    ss_X << X;
    std::cerr << "X:\n" << ss_X.str() << "\n";
    std::stringstream ss_Y;
    ss_Y << Y;
    std::cerr << "Y:\n" << ss_Y.str() << "\n";
#endif
        } else {
            // too far from a good fit, return NAs.
            stoat::LOG_WARN("Logistic regression didn't converge to a great fit (max score coeff " + 
                std::to_string(max_score) + ", max delta " + std::to_string(max_delta) +
                ") too far from a good fit. Skipping.", "count_number_not_convergence_skipping");
            return std::nan("");
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
        double XtWX_det = XtWX.determinant();
        if (XtWX_det <= .00000000001) {
            XtWX_det = .00000000001;
        }
        loglik += log(XtWX_det) / 2;
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
    Eigen::VectorXd beta0 = Eigen::VectorXd::Constant(X0.cols(), 0);
    Eigen::VectorXd delta0;
    Eigen::VectorXd score0;
    Eigen::MatrixXd XtWX0;
    Eigen::MatrixXd XtWXi0;

#ifdef DEBUG_STATS_TEST
    std::cerr << "X0:\n" << X0 << "\n";
    std::cerr << "beta:\n" << beta0 << "\n";
#endif

    iter = 0;
    while(iter < max_iterations) {
        Ypred = sigmoid(X0 * beta0);
        W = Eigen::MatrixXd::Constant(Y.size(), Y.size(), 0);
        for (size_t ii = 0; ii < Y.size(); ii++) {
            W(ii, ii) = Ypred(ii) * (1 - Ypred(ii));
        }
        XtWX0 = X0.transpose() * W * X0;
        XtWXi0 = inverse(XtWX0);
        if (penalize) {
            Eigen::MatrixXd sqrtW = Eigen::MatrixXd::Constant(Y.size(), Y.size(), 0);
            for (size_t ii = 0; ii < Y.size(); ii++) {
                sqrtW(ii, ii) = std::sqrt(Ypred(ii) * (1 - Ypred(ii)));
            }
            Eigen::MatrixXd hat = (sqrtW * X0) * XtWXi0 * (X0.transpose() * sqrtW);
            // the penalization uses the diagonal of this hat matrix as penalty
            // update the Ypred used for the score below with that penalty
            for (size_t ii = 0; ii < Y.size(); ii++) {
                Ypred(ii) = Ypred(ii) - (hat(ii, ii) * (0.5 - Ypred(ii)));
            }
        }
        score0 = X0.transpose() * (Y - Ypred);
        // first make sure the update is not too big (usually a bad sign)
        delta0 = XtWXi0 * score0;
        max_delta = delta0.cwiseAbs().maxCoeff();
        if (max_delta > cur_max_step) {
            // reduce the delta so that its largest component is cur_max_step
            delta0 = cur_max_step * delta0 / max_delta;
        }
        // update the betas
        beta0 = beta0 + delta0;
        iter++;

        // stop if it looks like we've converged
        if (score0.cwiseAbs().maxCoeff() < conv_tol && max_delta < conv_tol) {
            break;
        }
    }

#ifdef DEBUG_STATS_TEST
    std::cerr << "reduced model beta:\n" << beta0 << "\n\n";
#endif
    
    // compute the reduced model
    // compute log-likelihood for reduced model
    Eigen::VectorXd Xbeta0 = X0 * beta0;
    double loglik0 = 0;
    for (size_t ii = 0; ii < Y.size(); ii++) {
        loglik0 += Y(ii) * Xbeta0(ii) - std::log(1 + std::exp(Xbeta0(ii)));
    }
    if (penalize) {
        // if we're using the Firth penalized regression, the log-likelihood needs to
        // include the penalty, here half the determinant of the Fisher information matrix (Xt W X)
        Ypred = sigmoid(Xbeta0);
        W = Eigen::MatrixXd::Constant(Y.size(), Y.size(), 0);
        for (size_t ii = 0; ii < Y.size(); ii++) {
            W(ii, ii) = Ypred(ii) * (1 - Ypred(ii));
        }
        XtWX0 = X0.transpose() * W * X0;
        double XtWX_det = XtWX0.determinant();
        if (XtWX_det <= .00000000001) {
            XtWX_det = .00000000001;
        }
        loglik0 += log(XtWX_det) / 2;
    }
#ifdef DEBUG_STATS_TEST
    std::cout << "log-likelihood reduced model:" << loglik0 << "\n\n";
#endif
    
    // compute -2 * log-likelihood ratio of those models
    double loglik_ratio = -2 * (loglik0 - loglik);
    // should follow a chi2 distribution with df_full-df_reduced degrees of freedom
    size_t df = X.cols() - X0.cols();

#ifdef DEBUG_STATS_TEST
    std::cerr << "Chi2: " << loglik_ratio << " with df=" << df << "\n\n";
#endif

    if (loglik_ratio < 0) {
        stoat::LOG_WARN("Negative log_likelihood ratio, likely due to issues fitting the full/reduced models. Skipping.", "count_number_negative_log_likelihood");

#ifdef DEBUG_STATS_TEST
    std::stringstream ss_X;
    ss_X << X;
    std::cerr << "X:\n" << ss_X.str() << "\n";
    std::stringstream ss_Y;
    ss_Y << Y;
    std::cerr << "Y:\n" << ss_Y.str() << "\n";
#endif

        return std::nan("");
    }

    // compute p-value (double first)
    boost::math::chi_squared_distribution<double> dist(df);
    double p_value = boost::math::cdf(boost::math::complement(dist, loglik_ratio));

    return p_value;
}

// ------------------------ Chi2 test ------------------------

double FisherChi2::chi2_2x2(const size_t& a, const size_t& b, const size_t& c, const size_t& d) const {

    // Row and column sums
    double row1 = static_cast<double>(a + b);
    double row2 = static_cast<double>(c + d);
    double col1 = static_cast<double>(a + c);
    double col2 = static_cast<double>(b + d);
    double total = row1 + row2;

    // Early check: zero row or column
    if (row1 == 0 || row2 == 0 || col1 == 0 || col2 == 0) {
        stoat::LOG_WARN("Chi2 2x2: row or column sum is zero, returning NaN", "chi2_zero");
        return std::nan("");
    }

    // Expected counts
    double expected_a = row1 * col1 / total;
    double expected_b = row1 * col2 / total;
    double expected_c = row2 * col1 / total;
    double expected_d = row2 * col2 / total;

    // Check for zero, infinity, or NaN
    // matis: I don't really know what we expect to output here, at least I add a warning for user/log
    if (expected_a <= 0.0 || expected_b <= 0.0 || 
        expected_c <= 0.0 || expected_d <= 0.0 ||
        !std::isfinite(expected_a) || !std::isfinite(expected_b) || 
        !std::isfinite(expected_c) || !std::isfinite(expected_d)) {
        stoat::LOG_WARN("Chi2 2x2: expected counts out of valid range (zero or overflow), returning NaN", "chi2_inf");
        return std::nan("");
    }

    double chi2_stat = 0;
    chi2_stat += std::pow((double)a - expected_a, 2) / expected_a;
    chi2_stat += std::pow((double)b - expected_b, 2) / expected_b;
    chi2_stat += std::pow((double)c - expected_c, 2) / expected_c;
    chi2_stat += std::pow((double)d - expected_d, 2) / expected_d;

    // Degrees of freedom for 2×2 table
    size_t df = 1;

    // --- Compute p-value (double first) ---
    boost::math::chi_squared_distribution<double> dist(df);
    double p_value = boost::math::cdf(boost::math::complement(dist, chi2_stat));

    return p_value;
}

// Check if the observed matrix is valid (no zero rows/columns)
double FisherChi2::chi2_2xN(const std::vector<size_t>& g0, const std::vector<size_t>& g1) const {

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

    if (total == 0) {
        stoat::LOG_WARN("Total equal to 0. Skipping.", "count_number_total_zero");
        return std::nan("");
    }

    if (row_total_0 == 0 || row_total_1 == 0) {
        stoat::LOG_WARN("One of the rows has a total of 0. Skipping.", "count_number_row_zero");
        return std::nan("");
    }

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
    double p_value = boost::math::cdf(boost::math::complement(dist, chi2));

    return p_value;
}

// ------------------------ Fisher exact test ------------------------


std::pair<double, double> FisherChi2::fisher_chi2(const std::vector<size_t>& g0, const std::vector<size_t>& g1) const {
    
    double chi2_p_value = std::nan("");
    double fastfisher_p_value = std::nan("");
    
    // compute Fisher's exact or Chi-squared test p-value
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
double LinearRegression::linear_regression(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const size_t num_predictors) const {
#ifdef DEBUG_STATS_TEST
    std::cerr << "X:\n" << X << "\n";
    std::cerr << "Y:\n" << Y << "\n";
#endif
    size_t num_samples = Y.size();
    size_t num_covariates = X.cols() - num_predictors - 1;
    stoat::LOG_TRACE("Linear regression. " + std::to_string(num_samples) + " samples, " + std::to_string(num_predictors) + " predictors, " + std::to_string(num_covariates) + " covariates.");

    int num_params_full = 1 + num_predictors + num_covariates; // intercept + predictors + covariates
    int df_denominator = num_samples - num_params_full;

    // JEAN problem can arise if we have less samples than variables
    // maybe we should skip those tests when they happen? For now, warning the user and recommending increasing -I
    if (df_denominator <= 0) {
        stoat::LOG_WARN("Too few samples (" + std::to_string(num_samples) + ") compared to alleles+covariates (" + std::to_string(num_params_full) + ") in this snarl. Skipping. Note: increasing the minimum number of individuals with -I could help avoid those issues and get more robust associations in general.", "count_number_few_sample");

#ifdef DEBUG_STATS_TEST
    std::stringstream ss_X;
    ss_X << X;
    std::cerr << "X:\n" << ss_X.str() << "\n";
    std::stringstream ss_Y;
    ss_Y << Y;
    std::cerr << "Y:\n" << ss_Y.str() << "\n";
#endif

        return std::nan("");
    }

    // fit the full model
    // Solve beta X = Y using QR decomposition
    // colPivHouseholderQr is faster than fullPivHouseholderQr but less stable. could switch to full if we encounter problem (or try the HouseolderQR which is less stable but even faster)
    Eigen::VectorXd beta = X.colPivHouseholderQr().solve(Y);
    Eigen::VectorXd y_hat = X * beta;

#ifdef DEBUG_STATS_TEST
    std::cerr << "beta:\n" << beta << "\n";
#endif

    // sum of squared errors and total sum of squares total (like using the mean as prediction)
    double SSE_full = 0.0;
    double SST = 0.0;
    double y_mean = Y.sum() / Y.size();
    for (int i = 0; i < num_samples; ++i) {
        SSE_full += (Y(i) - y_hat(i)) * (Y(i) - y_hat(i));
        SST += (Y(i) - y_mean) * (Y(i) - y_mean);
    }

    // reduced model without any predictor of interest (allele counts)
    double SSE_reduced = 0.0;
    if (num_covariates > 0) {
        // keep the intercept and covariates only
        Eigen::MatrixXd X_reduced = Eigen::MatrixXd::Constant(X.rows(), 1 + num_covariates, 1);
        X_reduced.block(0, 1, X.rows(), num_covariates) = X.block(0, 1 + num_predictors, X.rows(), num_covariates);

        #ifdef DEBUG_STATS_TEST
        std::cerr << "X reduced:\n" << X_reduced << "\n";
#endif

        Eigen::VectorXd beta_r = X_reduced.colPivHouseholderQr().solve(Y);
        Eigen::VectorXd y_hat_r = X_reduced * beta_r;
        for (int i = 0; i < num_samples; ++i) {
            SSE_reduced += (Y(i) - y_hat_r(i)) * (Y(i) - y_hat_r(i));
        }
    } else {
        // no covariates: reduced model = intercept only
        SSE_reduced = SST;
    }

    if (SSE_full == 0) {
        // rare but can happen. What should we do?
        if (SSE_reduced == 0) {
            stoat::LOG_WARN("SSE is null for both the full and reduced model. Skipping.", "count_number_sse_null");

#ifdef DEBUG_STATS_TEST
    std::stringstream ss_X;
    ss_X << X;
    std::cerr << "X:\n" << ss_X.str() << "\n";
    std::stringstream ss_Y;
    ss_Y << Y;
    std::cerr << "Y:\n" << ss_Y.str() << "\n";
#endif

            return std::nan("");
        } else {
            // not sure what to do here. Adding 0.1% errors to the perfect predictions?
            stoat::LOG_DEBUG("SSE_full == 0. Adding 0.1 percentage errors to the perfect predictions.");
            for (int i = 0; i < num_samples; ++i) {
                SSE_full += 0.001 * 0.001;
            }
        }
    }
    
    // F = [(SSE_reduced - SSE_full) / df_numerator] / [SSE_full / df_denominator]
    // Numerator df = number of tested predictors
    size_t df_numerator = num_predictors;
    double numerator = (SSE_reduced - SSE_full) / df_numerator;

    // Denominator df = residual df in full model
    double denominator = SSE_full / df_denominator;

    // time to compute the F-statistic
    double F_stat = numerator / denominator;

    if (F_stat < 0 && F_stat > -0.00001){
        stoat::LOG_WARN("F statistic is negative but very close to 0 (" + std::to_string(F_stat) + "). Assuming it's 0.", "count_number_f_stat_close_to_0");

#ifdef DEBUG_STATS_TEST
    std::stringstream ss_X;
    ss_X << X;
    std::cerr << "X:\n" << ss_X.str() << "\n";
    std::stringstream ss_Y;
    ss_Y << Y;
    std::cerr << "Y:\n" << ss_Y.str() << "\n";
#endif

        F_stat = 0;
    }

    if (F_stat < -0.00001){
        stoat::LOG_WARN("F statistic is negative: " + std::to_string(F_stat) + " = " + std::to_string(numerator) + "/" + std::to_string(denominator) + ". This is concerning, skipping. Recommendation: increase the minimum number of individuals with -I to get more robust associations and avoid issues.", "count_number_f_stat_negative");
        stoat::LOG_DEBUG("SSE_reduced: " + std::to_string(SSE_reduced) + ",SSE_full: " + std::to_string(SSE_full));

#ifdef DEBUG_STATS_TEST
    std::stringstream ss_X;
    ss_X << X;
    std::cerr << "X:\n" << ss_X.str() << "\n";
    std::stringstream ss_Y;
    ss_Y << Y;
    std::cerr << "Y:\n" << ss_Y.str() << "\n";
#endif

        return std::nan("");
    }

    // compute a P-value
    boost::math::fisher_f_distribution<double> dist(df_numerator, df_denominator);
    double p_value = boost::math::cdf(boost::math::complement(dist, F_stat));

    return p_value;
}

} // namespace stoat
