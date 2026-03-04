#include <Eigen/Dense>
#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>
#include <vector>
#include <cmath>

Eigen::MatrixXd pseudoInverse(const Eigen::MatrixXd& A, double tol = 1e-10) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        A, Eigen::ComputeThinU | Eigen::ComputeThinV
    );

    const auto& U = svd.matrixU();
    const auto& V = svd.matrixV();
    const auto& S = svd.singularValues();

    Eigen::MatrixXd S_pinv = Eigen::MatrixXd::Zero(V.cols(), U.cols());

    for (int i = 0; i < S.size(); ++i) {
        if (S(i) > tol)
            S_pinv(i, i) = 1.0 / S(i);
    }

    return V * S_pinv * U.transpose();
}

// Try normal inverse if well-conditioned, otherwise use pseudo-inverse
Eigen::MatrixXd invertOrPseudoInvert(const Eigen::MatrixXd& A, double tol = 1e-10) {
    // Sanity check
    if (A.rows() != A.cols()) {
        // Non-square → pseudo-inverse only
        std::cout << "Warning: Matrix is not square (rows=" << A.rows() << ", cols=" << A.cols() << "). Using pseudo-inverse.\n";
        return pseudoInverse(A, tol);
    }

    // Try LDLT decomposition (fast & stable for symmetric matrices)
    Eigen::LDLT<Eigen::MatrixXd> ldlt(A);

    // Check if matrix is numerically invertible
    if (ldlt.info() == Eigen::Success) {
        // Check conditioning via diagonal D
        const auto& D = ldlt.vectorD();
        for (int i = 0; i < D.size(); ++i) {
            if (std::abs(D(i)) < tol) {
                // Near-singular
                std::cout << "Warning: Matrix is near-singular (D[" << i << "]=" << D(i) << "). Using pseudo-inverse.\n";
                return pseudoInverse(A, tol);
            }
        }

        // Simple inversion
        return ldlt.solve(Eigen::MatrixXd::Identity(A.rows(), A.cols()));
    }

    // Fallback
    std::cout << "Warning: Matrix is not invertible. Using pseudo-inverse.\n";
    return pseudoInverse(A, tol);
}

// Log-likelihood for logistic regression
double log_likelihood(const Eigen::MatrixXd& X,
                      const Eigen::VectorXd& y,
                      const Eigen::VectorXd& beta) {

    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = 1.0 / (1.0 + (-eta.array()).exp());

    const double eps = 1e-12; // avoid log(0)
    return (y.array() * (mu.array() + eps).log() +
            (1 - y.array()) * (1 - mu.array() + eps).log()).sum();
}

// IRLS solver
// Logistic regression using pure Newton-Raphson
Eigen::VectorXd fit_logistic_newton(const Eigen::MatrixXd& X,
                                    const Eigen::VectorXd& y,
                                    int max_iter = 25,
                                    double tol = 0.001) {

    Eigen::VectorXd beta = Eigen::VectorXd::Zero(X.cols());
    Eigen::VectorXd eta, mu, w, gradient, delta;
    Eigen::MatrixXd Xt = X.transpose();
    Eigen::MatrixXd XtWX;

    for (int iter = 0; iter < max_iter; ++iter) {

        std::cout << "Iteration " << iter << "\n";

        // η = Xβ
        eta = X * beta;
        std::cout << "Eta = " << eta.transpose() << "\n";

        // μ = sigmoid(η)
        mu = (1.0 + (-eta.array()).exp()).inverse();
        std::cout << "Mu = " << mu.transpose() << "\n";

        // W = diag( μ(1-μ) ), here it's only μ(1-μ) them we will use as a diagonal matrix
        w = mu.array() * (1.0 - mu.array());
        std::cout << "W = " << w.transpose() << "\n";

        // XtWX = Xᵀ W X
        XtWX = Xt * w.asDiagonal() * X;
        std::cout << "XtWX = \n" << XtWX << "\n";

        // Gradient: Xᵀ(y - μ)
        gradient = Xt * (y - mu);
        std::cout << "Gradient = " << gradient.transpose() << "\n";

        // Solve: (Xᵀ W X) Δβ = gradient
        delta = invertOrPseudoInvert(XtWX) * gradient;
        std::cout << "Delta = " << delta.transpose() << "\n";

        // Update
        beta += delta;
        std::cout << "Beta = " << beta.transpose() << "\n";

        // Convergence check: ||Δβ|| < tol
        if (delta.norm() < tol)
            break;
    }

    return beta;
}

// Logistic regression with Likelihood Ratio Test
void logistic_regression(const Eigen::MatrixXd& X,
                         const Eigen::VectorXd& y,
                         const size_t num_predictors,
                         const size_t num_covariates) {

    // -----------------------------
    // Full model
    // -----------------------------
    Eigen::VectorXd beta_full = fit_logistic_newton(X, y);
    double ll_full = log_likelihood(X, y, beta_full);

    // -----------------------------
    // Reduced model (intercept + covariates only)
    // -----------------------------
    Eigen::MatrixXd X_reduced(X.rows(), 1 + num_covariates);

    // 1. Intercept
    X_reduced.col(0) = X.col(0);

    // 2. Covariates (they start after predictors)
    X_reduced.rightCols(num_covariates) = X.middleCols(1 + num_predictors, num_covariates);
    Eigen::VectorXd beta_reduced = fit_logistic_newton(X_reduced, y);
    double ll_reduced = log_likelihood(X_reduced, y, beta_reduced);

    // -----------------------------
    // Likelihood Ratio Test
    // -----------------------------
    double LRT_stat = 2.0 * (ll_full - ll_reduced);

    if (LRT_stat < 0) {
        std::cerr << "Warning: LRT statistic is negative (ll_full < ll_reduced). Setting LRT_stat to 0.\n";
        LRT_stat = 0; // numerical safeguard
    }

    // Degrees of freedom = difference in number of parameters
    int df = X.cols() - X_reduced.cols();

    boost::math::chi_squared chi_sq(df);
    double p_value = 1.0 - boost::math::cdf(chi_sq, LRT_stat);

    // -----------------------------
    // Print results
    // -----------------------------
    std::cout << "Logistic Regression Likelihood Ratio Test\n";
    std::cout << "------------------------------------------\n";
    std::cout << "Log-likelihood (Full):    " << ll_full << "\n";
    std::cout << "Log-likelihood (Reduced): " << ll_reduced << "\n";
    std::cout << "LRT statistic:            " << LRT_stat << "\n";
    std::cout << "Degrees of freedom:       " << df << "\n";
    std::cout << "p-value:                  " << p_value << "\n";
}

int main() {

    Eigen::MatrixXd X(10, 7);
    X <<
        1, 0, 0.5, 0, 10,  2,  5,
        1, 1, 0, 0, 8,  1,  3,
        1, 0, 1, 0,  9,  4,  2,
        1, 1, 0, 0,  7,  3,  1,
        1, 0, 0.5, 0,  11,  2,  4,
        1, 1, 0, 0,  6,  1,  2,
        1, 0, 1, 0, 10,  3,  5,
        1, 1, 0, 0, 9,  2,  1,
        1, 0, 0, 0, 8,  4,  3,
        1, 0.5, 0, 0.5, 7,  2,  4;

    Eigen::VectorXd y(10);
    y << 1, 1, 1, 1, 1, 0, 0, 0, 0, 0;

    logistic_regression(X, y, 3, 3);

    return 0;
}

// g++ -O2 -std=c++17 logistic_regression_matis.cpp -I /usr/include/eigen3 -o logistic_regression_matis && ./logistic_regression_matis
