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

// if the matrix looks non-invertible, use the pseudo-inverse (see above)
Eigen::MatrixXd inverse(const Eigen::MatrixXd &A) {
    if (A.determinant() < 1e-6) {
        std::cout << "Warning: Matrix is near-singular (determinant=" << A.determinant() << "). Using pseudo-inverse.\n";
        return pseudoInverse(A);
    } else {
        return A.inverse();
    }
}

// compute the sigmoid function, here corresponding to the "predicted" probability associated to a set of observations and the model (Beta x X)
Eigen::VectorXd sigmoid(const Eigen::VectorXd& t) {
    Eigen::VectorXd res = Eigen::VectorXd::Constant(t.size(), 0);
    for (size_t idx = 0; idx < t.size(); idx++) {
        res(idx) = 1 / (1 + std::exp(-t(idx)));
    }
    return res;
}

// logistic regression using the Maximum Likelihood Estimate with Newton-Raphson method
// Logistic regression with Likelihood Ratio Test
void logistic_regression(const Eigen::MatrixXd& X,
                         const Eigen::VectorXd& Y,
                         const size_t num_predictors,
                         const size_t num_covariates) {

    // prepare an output objet and init to NA
    std::string pv = "NA";
    size_t max_iterations = 25;
    double conv_tol = 0.001;
    double max_step = 3;

    size_t num_samples = Y.size();
    size_t num_features = X.cols();

    // we'll look for the beta coefficients that maximize the likelihood using the Newton-Raphson method
    // it's an iterative approach so we start by initializing the coefficients (to 0)
    Eigen::VectorXd beta = Eigen::VectorXd::Constant(X.cols(), 0);

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
        std::cout << "Iteration " << iter << "\n";
        // compute the predicted probs from the current beta coefficients
        Ypred = sigmoid(X * beta);

        std::cout << "Ypred = " << Ypred.transpose() << "\n";
        // W is a diag matrix with S(BXi)(1-S(BXi)) on the diag
        // JEAN might be quite big and used for the diag. Could do better here (sparse matrix? manually multiplying X's columns?)
        W = Eigen::MatrixXd::Constant(Y.size(), Y.size(), 0);
        for (size_t ii = 0; ii < Y.size(); ii++) {
            W(ii, ii) = Ypred(ii) * (1 - Ypred(ii));
        }

        std::cout << "W = " << W.transpose() << "\n";
    
        XtWX = Xt * W * X;
        XtWXi = inverse(XtWX);
        // if we want to penalize compute the penalty for the current betas
        if (penalize) {
            std::cout << "Penalizing iteration " << iter << " with max step " << cur_max_step << "\n";
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
        std::cout << "Score = " << score.transpose() << "\n";

        // update the  beta coefficients using the Newton-Raphton formula (XtWX is the hessian)
        // first make sure the update is not too big (usually a bad sign)
        delta = XtWXi * score;
        std::cout << "Delta = " << delta.transpose() << "\n";

        max_delta = delta.cwiseAbs().maxCoeff();
        if (max_delta > cur_max_step) {
            std::cout << "Max delta " << max_delta << " is greater than current max step " << cur_max_step << ". Reducing delta.\n";
            // reduce the delta so that its largest component is cur_max_step
            delta = cur_max_step * delta / max_delta;
            std::cout << "Reduced delta = " << delta.transpose() << "\n";
        }
        // ready to update the betas
        beta += delta;
        std::cout << "Beta = " << beta.transpose() << "\n";

        // switch to penalizing if betas get too large (and we were not already penalizing)
        // JEAN any beta larger than 10 is suspicious and unexpected, but maybe there is a better way to decide if we should switch to Firth's regression
        if (!penalize && (beta.maxCoeff() > 10 || beta.maxCoeff() < -10)) {
            std::cout << "Penalization activation" << std::endl;
            penalize = true;
            // start again the iteration process
            iter = 0;
            beta = Eigen::VectorXd::Constant(X.cols(), 0);

            continue;
        } else {
            iter++;
        }

        max_score = score.cwiseAbs().maxCoeff();
        std::cout << "Max delta = " << max_delta << "\n";
        std::cout << "Max score = " << max_score << "\n";
        std::cout << "conv_tol = " << conv_tol << "\n";
        // stop if it looks like we've converged (small score and small deltas)
        if (max_score < conv_tol && max_delta < conv_tol) {
            std::cout << "Convergence achieved with max score " << max_score << " and max delta " << max_delta << ". Stopping iterations.\n";
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
            std::cout << "Logistic regression didn't converge to a great fit (max score coeff " + std::to_string(max_score) + ", max delta " + std::to_string(max_delta) + ") but not too bad. Continuing.";
        } else {
            // too far from a good fit, return NAs.
            std::cout << "Logistic regression didn't converge to a good fit (max score coeff " + std::to_string(max_score) + ", max delta " + std::to_string(max_delta) + "). Returning NA.";
            return;
        }
    }

    // compute log-likelihood for current betas
    Eigen::VectorXd Xbeta = X * beta;
    double loglik = 0;
    for (size_t ii = 0; ii < Y.size(); ii++) {
        loglik += Y(ii) * Xbeta(ii) - std::log(1 + std::exp(Xbeta(ii)));
    }
    if (penalize) {
        std::cout << "Penalizing n°2 iteration " << iter << " with max step " << cur_max_step << "\n";
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
        beta0 += delta0;
        iter++;

        // stop if it looks like we've converged
        if (score0.cwiseAbs().maxCoeff() < conv_tol && max_delta < conv_tol) {
            break;
        }
    }
    
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
    
    // compute -2 * log-likelihood ratio of those models
    double loglik_ratio = -2 * (loglik0 - loglik);
    // should follow a chi2 distribution with df_full-df_reduced degrees of freedom
    int df = X.cols() - X0.cols();

    if (loglik_ratio < 0) {
        std::cout <<"Negative log_likelihood ratio, likely due to issues fitting the full/reduced models. Skipping.";
        return;
    }
    
    // compute p-value (double first)
    boost::math::chi_squared_distribution<double> dist(df);
    double p_value = 1.0 - boost::math::cdf(dist, loglik_ratio);

    // -----------------------------
    // Print results
    // -----------------------------
    std::cout << "Logistic Regression Likelihood Ratio Test\n";
    std::cout << "------------------------------------------\n";
    std::cout << "Log-likelihood (Full):    " << loglik << "\n";
    std::cout << "Log-likelihood (Reduced): " << loglik0 << "\n";
    std::cout << "loglik_ratio:             " << loglik_ratio << "\n";
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

// g++ -O2 -std=c++17 logistic_regression_jean.cpp -I /usr/include/eigen3 -o logistic_regression_jean && ./logistic_regression_jean
