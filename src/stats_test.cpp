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

// ------------------------ Filtration Function ------------------------

// Return true when snarl must be filtered and false if not
bool filtration_quantitative_table(
    const std::vector<std::vector<double>>& X,
    const size_t& min_individuals,
    const double& maf_threshold) {
    
    // number of path < 2 OR not enougth individuals
    if (X.empty() || X[0].size() < 2 || X.size() < min_individuals) {
        if (X.empty()) {
            stoat::LOG_DEBUG("Filtration cause: X is empty.");
            return true;
        }

        if (X[0].size() < 2) {
            stoat::LOG_DEBUG("Filtration cause: Not enough paths (" + std::to_string(X[0].size()) + " < 2)");
            return true;
        }

        if (X.size() < min_individuals) {
            stoat::LOG_DEBUG("Filtration cause: Not enough individuals (" + std::to_string(X.size()) + " < " + std::to_string(min_individuals) + ")");
            return true;
        }
    }

    size_t numPaths = X[0].size();
    std::vector<double> table(numPaths, 0.0);
    double totalSum = 0.0;

    // Compute column sums and total sum
    for (const auto& row : X) {
        for (size_t i = 0; i < numPaths; ++i) {
            table[i] += row[i];
            totalSum += row[i];
        }
    }

    int count_above_threshold = 0;
    for (size_t i = 0; i < numPaths; ++i) {
        double freq = table[i] / totalSum;
        double maf = std::min(freq, 1.0 - freq);
        if (maf > maf_threshold) {
            ++count_above_threshold;
        }
    }

    return count_above_threshold < 2;
}

void remove_empty_columns_quantitative_table(
    std::vector<std::vector<double>>& X) {

    if (X.empty() || X[0].size() < 2) return;

    size_t num_rows = X.size();
    size_t num_cols = X[0].size();

    // Identify non-empty columns
    std::vector<bool> keep_column(num_cols, false);

    for (size_t col = 0; col < num_cols; ++col) {
        for (size_t row = 0; row < num_rows; ++row) {
            double val = X[row][col];
            if (val != 0.0 && !std::isnan(val)) {
                keep_column[col] = true;
                break;
            }
        }
    }

    // Create filtered X
    std::vector<std::vector<double>> df_filtered;
    df_filtered.reserve(num_rows);

    for (size_t row = 0; row < num_rows; ++row) {
        std::vector<double> new_row;
        for (size_t col = 0; col < num_cols; ++col) {
            if (keep_column[col]) {
                new_row.push_back(X[row][col]);
            }
        }
        df_filtered.push_back(std::move(new_row));
    }

    // Replace original X with filtered one
    X = std::move(df_filtered);
}

bool check_last_columns_quantitative_table(
    const std::vector<std::vector<double>>& X) {

    if (X[0].size() > 1) return false;

    size_t num_rows = X.size();

    // Check if the lonely columns have identical values
    for (size_t r = 1; r < num_rows-1; ++r) {
        if (X[r][0] != X[0][0]) {
            return false; // Not identical, keep
        }
    }

    return true; // Identical, filter out
}

void combine_identical_columns_quantitative_table(
    std::vector<std::vector<double>>& X) {

    size_t num_rows = X.size();
    size_t num_cols = X[0].size();

    if (num_cols < 3) {return;} // avoid creation of unique column

    std::vector<bool> merged(num_cols, false);
    std::vector<std::vector<double>> new_cols;

    for (size_t i = 0; i < num_cols; ++i) {
        if (merged[i]) continue;

        std::vector<double> new_col = X[0][i] == X[0][i] ? std::vector<double>(num_rows, 0.0) : X[0]; // Init new column

        // Start with current column
        for (size_t r = 0; r < num_rows; ++r) {
            new_col[r] = X[r][i];
        }

        // Try to find identical columns
        for (size_t j = i + 1; j < num_cols; ++j) {
            if (merged[j]) continue;

            bool identical = true;
            for (size_t r = 0; r < num_rows; ++r) {
                if (X[r][j] != X[r][i]) {
                    identical = false;
                    break;
                }
            }

            // If identical, sum the column into new_col
            if (identical) {
                for (size_t r = 0; r < num_rows; ++r) {
                    new_col[r] += X[r][j];
                }
                merged[j] = true;
            }
        }

        new_cols.push_back(std::move(new_col));
    }

    // Rebuild X from new_cols (transpose)
    std::vector<std::vector<double>> result(num_rows, std::vector<double>(new_cols.size()));
    for (size_t r = 0; r < num_rows; ++r) {
        for (size_t c = 0; c < new_cols.size(); ++c) {
            result[r][c] = new_cols[c][r];
        }
    }

    X = std::move(result);
}

void remove_last_columns_quantitative_table(std::vector<std::vector<double>>& X) {
    if (X.empty() || X[0].empty()) return;

    for (auto& row : X) {
        if (!row.empty()) {
            row.pop_back(); // Remove last column from each row
        }
    }
}

void remove_empty_columns_binary_table(
    std::vector<size_t>& g0, 
    std::vector<size_t>& g1) {

    std::vector<size_t> g0_filtered;
    std::vector<size_t> g1_filtered;

    for (size_t i = 0; i < g0.size(); ++i) {
        if (g0[i] + g1[i] != 0) {
            g0_filtered.push_back(g0[i]);
            g1_filtered.push_back(g1[i]);
        }
    }

    g0 = std::move(g0_filtered);
    g1 = std::move(g1_filtered);
}

// true : filtration on; false : no filtration
bool filtration_binary_table(
    std::vector<size_t>& g0, 
    std::vector<size_t>& g1,
    const size_t& individuals_included, 
    const size_t& min_individuals,
    const double& maf_threshold) {

    // Not enougth individuals OR not enougth haplotypes OR number of paths < 2
    if (individuals_included < min_individuals || g0.size() < 2) {
        return true; // Empty or invalid input → filter
    }

    size_t haplotype_count = 0;
    for (size_t i = 0 ; i < g0.size() ; i++) {
        haplotype_count += g0[i];
        haplotype_count += g1[i];
    }

    int count_above_threshold = 0;
    for (size_t i = 0; i < g0.size(); ++i) {
        size_t columnSum = g0[i] + g1[i];

        double freq1 = static_cast<double>(columnSum) / haplotype_count;
        double maf = std::min(freq1, 1.0 - freq1);

        if (maf > maf_threshold) {
            ++count_above_threshold;
        }
    }

    return count_above_threshold < 2; // Keep if at least two MAFs path > MAF threshold
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
std::tuple<std::string, std::string, std::string> LogisticRegression::logistic_regression(
    const std::vector<std::vector<double>>& df,
    const std::vector<bool>& phenotype,
    const std::vector<std::vector<double>>& covariates) {

    size_t num_samples = df.size();
    size_t num_variants = df[0].size();
    size_t num_covariates = 0;
    size_t num_features = num_variants + 1; // +1 for intercept

    if (!covariates.empty()) {
        num_covariates = covariates[0].size();
        num_features = num_variants + num_covariates + 1; // +1 for intercept
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

// Multiply matrix A (m×n) with matrix B (n×p) -> result is m×p
std::vector<std::vector<double>> LinearRegression::matmul(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B) {
    int m = A.size(), n = A[0].size(), p = B[0].size();
    std::vector<std::vector<double>> result(m, std::vector<double>(p, 0.0));
    for (int i = 0; i < m; ++i)
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < p; ++j)
                result[i][j] += A[i][k] * B[k][j];
    return result;
}

// Multiply matrix A (m×n) with vector b (n) -> result is vector of size m
std::vector<double> LinearRegression::matvec(const std::vector<std::vector<double>> &A, const std::vector<double> &b) {
    int m = A.size(), n = A[0].size();
    std::vector<double> result(m, 0.0);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            result[i] += A[i][j] * b[j];
    return result;
}

// Transpose of a matrix
std::vector<std::vector<double>> LinearRegression::transpose(const std::vector<std::vector<double>> &A) {
    int m = A.size(), n = A[0].size();
    std::vector<std::vector<double>> result(n, std::vector<double>(m, 0.0));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            result[j][i] = A[i][j];
    return result;
}

// Convert std::vector<std::vector<double>> to Eigen::MatrixXd
Eigen::MatrixXd LinearRegression::toEigenMatrix(const std::vector<std::vector<double>>& mat) {
    int rows = mat.size();
    int cols = mat[0].size();
    Eigen::MatrixXd result(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result(i, j) = mat[i][j];
    return result;
}

// Convert Eigen::MatrixXd back to std::vector<std::vector<double>>
std::vector<std::vector<double>> LinearRegression::fromEigenMatrix(const Eigen::MatrixXd& mat) {
    int rows = mat.rows();
    int cols = mat.cols();
    std::vector<std::vector<double>> result(rows, std::vector<double>(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result[i][j] = mat(i, j);
    return result;
}

// Compute Moore-Penrose pseudoinverse using SVD
std::vector<std::vector<double>> LinearRegression::pseudoInverse(const std::vector<std::vector<double>>& A, double tol) {
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
std::vector<std::vector<double>> LinearRegression::inverse(
    const std::vector<std::vector<double>> &A, double tol ) {

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
            LOG_DEBUG("Using pseudo-inverse.");
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

// Performs linear regression and F-test for predictors only
std::tuple<std::string, std::string> LinearRegression::linear_regression(
    const std::vector<std::vector<double>>& X_raw,
    const std::vector<double>& y,
    const std::vector<std::vector<double>>& covariates) {

    int num_samples = X_raw.size();                                     // number of observations
    int num_predictors = X_raw[0].size();                               // number of predictors
    int num_covariates = covariates.empty() ? 0 : covariates[0].size(); // number of covariates

    LOG_TRACE("Linear regression. " + std::to_string(num_samples) + " samples, " + std::to_string(num_predictors) + " predictors, " + std::to_string(num_covariates) + " covariates.");
    
    // ---- FULL MODEL ----
    int num_params_full = 1 + num_predictors + num_covariates; // intercept + predictors + covariates
    std::vector<std::vector<double>> X_full(num_samples, std::vector<double>(num_params_full, 1.0));

    for (int i = 0; i < num_samples; ++i) {
        // Copy predictor variables
        std::copy(X_raw[i].begin(), X_raw[i].end(), X_full[i].begin() + 1);

        // Copy covariates if present
        if (num_covariates > 0) {
            std::copy(covariates[i].begin(), covariates[i].end(), X_full[i].begin() + 1 + num_predictors);
        }
    }

    // ---- Fit FULL model ----
    auto Xt = transpose(X_full);
    auto XtX = matmul(Xt, X_full);
    auto XtXi = inverse(XtX);
    auto Xty = matvec(Xt, y);
    auto beta = matvec(XtXi, Xty);
    auto y_hat = matvec(X_full, beta);

    // SSE (Sum of Squared Errors) for full model
    double SSE_full = 0.0;
    for (int i = 0; i < num_samples; ++i) {SSE_full += (y[i] - y_hat[i]) * (y[i] - y_hat[i]);}

    // ---- Compute R² ----
    double y_mean = std::accumulate(y.begin(), y.end(), 0.0) / num_samples;
    double SST = 0.0;
    for (int i = 0; i < num_samples; ++i) {SST += (y[i] - y_mean) * (y[i] - y_mean);}
    double R2 = 1.0 - SSE_full / SST;

    // ---- REDUCED MODEL ----
    double SSE_reduced = 0.0;
    if (num_covariates > 0) {
        // reduced model: intercept + covariates only
        int num_params_reduced = 1 + num_covariates;
        std::vector<std::vector<double>> X_reduced(num_samples, std::vector<double>(num_params_reduced, 1.0));

        for (int i = 0; i < num_samples; ++i) {
            for (int j = 0; j < num_covariates; ++j) {
                X_reduced[i][1 + j] = covariates[i][j];
            }
        }

        auto Xt_r = transpose(X_reduced);
        auto XtX_r = matmul(Xt_r, X_reduced);
        auto XtXi_r = inverse(XtX_r);
        auto Xty_r = matvec(Xt_r, y);
        auto beta_r = matvec(XtXi_r, Xty_r);
        auto y_hat_r = matvec(X_reduced, beta_r);

        for (int i = 0; i < num_samples; ++i) {SSE_reduced += (y[i] - y_hat_r[i]) * (y[i] - y_hat_r[i]);}
    } else {
        // no covariates: reduced model = intercept only
        // JEAN same as SST already computed above then, right?
        // for (int i = 0; i < num_samples; ++i) {SSE_reduced += (y[i] - y_mean) * (y[i] - y_mean);}
        SSE_reduced = SST;
    }

    // ---- F-statistic ----
    // Numerator df = number of tested predictors
    // Denominator df = residual df in full model
    int df_numerator = num_predictors;
    int df_denominator = (num_samples - num_params_full) <= 0 ? 1 : num_samples - num_params_full;

    // Compute F-statistic:
    // F = [(SSE_reduced - SSE_full) / df_numerator] / [SSE_full / df_denominator]
    double numerator = (SSE_reduced - SSE_full) / df_numerator;
    double denominator = SSE_full / df_denominator;
    double F_stat = numerator / denominator;

    // JEAN problem can arise if we have less samples than variables
    // maybe we should skip those tests when they happen? For now, warning the user and recommending increasing -I
    if (num_samples < num_params_full) {
        stoat::LOG_WARN("Less samples than alleles+covariates in some snarls. Increasing the minimum number of individuals with -I to get more robust associations and avoid issues.");
    }
    assert(F_stat > 0);
    boost::math::fisher_f_distribution<double> dist(df_numerator, df_denominator);
    double p_value = 1.0 - boost::math::cdf(dist, F_stat);

    std::string p_value_str;

    if (p_value == 0.0 || !std::isfinite(p_value)) {

        // Recompute only if underflow or invalid result
        const boost::multiprecision::cpp_dec_float_50 F_hp      = F_stat;
        const boost::multiprecision::cpp_dec_float_50 df_n_hp   = df_numerator;
        const boost::multiprecision::cpp_dec_float_50 df_d_hp   = df_denominator;

        boost::math::fisher_f_distribution<boost::multiprecision::cpp_dec_float_50> dist_hp(df_n_hp, df_d_hp);
        boost::multiprecision::cpp_dec_float_50 p_hp = boost::multiprecision::cpp_dec_float_50(1) - boost::math::cdf(dist_hp, F_hp);

        p_value_str = stoat::set_precision_float_50(p_hp);
    
    } else {
        p_value_str = stoat::set_precision(p_value);
    }

    // Format R² value
    const std::string r2_str = stoat::set_precision(R2);
    return std::make_tuple(p_value_str, r2_str);
}

} // namespace stoat
