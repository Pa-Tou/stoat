#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using boost::multiprecision::cpp_dec_float_50;
using boost::math::chi_squared_distribution;

// In the issue #53 : https://github.com/Pa-Tou/stoat/issues/53
// We discovert that collinearity still could create problem over linear regression
// computing, so my idea is to try to do like R does : Identify depend column and 
// remove them from the X, Here we should conserve the predicator and the intercept column
// to evaluate depend column I used Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X) & qr.setThreshold(1e-7);
int main() {

    Eigen::MatrixXd X(8,7);

    // INTERCEPT + PREDICATOR MUST STAY
    std::vector<size_t> predicator_index = {0,1};
    size_t num_predictors = 1;

    std::unordered_set<size_t> protected_cols(
        predicator_index.begin(),
        predicator_index.end()
    );

    X << 1,0,0,0,0,0,1,
         1,1,1,1,0,0,2,
         1,0,0,0,0,0,1,
         1,0,1,1,0,0,1,
         1,1,0,0,1,0,2,
         1,0,0,0,0,1,2,
         1,0,1,0,0,0,1,
         1,0,0,0,0,0,1;

    Eigen::VectorXd Y(8);

    Y << 39.347,
        53.6,
        60.8,
        54.596,
        70.351,
        51.429,
        64.16,
        44.082;

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
    qr.setThreshold(1e-7);

    int rank = qr.rank();
    Eigen::VectorXi perm = qr.colsPermutation().indices();

    std::vector<int> keep;
    std::vector<int> drop;

    for(int i = 0; i < X.cols(); i++) {
        if(i < rank)
            keep.push_back(perm(i));
        else
            drop.push_back(perm(i));
    }

    // enforce protected columns
    for(size_t p : predicator_index) {
        if(std::find(drop.begin(), drop.end(), p) != drop.end()) {
            for(size_t i = 0; i < keep.size(); i++) {
                if(!protected_cols.count(keep[i])) {
                    drop.erase(std::remove(drop.begin(), drop.end(), p), drop.end());
                    drop.push_back(keep[i]);
                    keep[i] = p;
                    break;
                }
            }
        }
    }

    std::sort(keep.begin(), keep.end());

    std::cout << "Final kept columns: ";
    for(auto c : keep)
        std::cout << c << " ";
    std::cout << std::endl;

    Eigen::MatrixXd X_filtered(X.rows(), keep.size());

    for(size_t j = 0; j < keep.size(); j++)
        X_filtered.col(j) = X.col(keep[j]);

    std::cout << "\nMatrix:\n" << X_filtered << std::endl;

    size_t num_samples = Y.size();
    size_t num_covariates = X_filtered.cols() - num_predictors - 1;

    std::cout << "Linear regression. " + std::to_string(num_samples) + " samples, " + std::to_string(num_predictors) + " predictors, " + std::to_string(num_covariates) + " covariates." << std::endl;

    // fit the full model
    // Solve beta X = Y using QR decomposition
    // colPivHouseholderQr is faster than fullPivHouseholderQr but less stable. could switch to full if we encounter problem (or try the HouseolderQR which is less stable but even faster)
    Eigen::VectorXd beta = X_filtered.colPivHouseholderQr().solve(Y);
    Eigen::VectorXd y_hat = X_filtered * beta;

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
        Eigen::MatrixXd X_reduced = Eigen::MatrixXd::Constant(X_filtered.rows(), 1 + num_covariates, 1);
        X_reduced.block(0, 1, X_filtered.rows(), num_covariates) = X_filtered.block(0, 1 + num_predictors, X_reduced.rows(), num_covariates);

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
            std::cout << "SSE is null for both the full and reduced model. Skipping."<< std::endl;
            std::stringstream ss;
            ss << "X_filtered:\n";
            for (int i = 0; i < X_filtered.rows(); ++i) {
                for (int j = 0; j < X_filtered.cols(); ++j) {
                    ss << X_filtered(i,j) << " ";
                }
                ss << "\n";
            }
            std::cout << ss.str()<< std::endl;

            std::stringstream ss_Y;
            ss_Y << "Y:\n";
            for (int i = 0; i < Y.size(); ++i) {
                ss_Y << Y[i] << " ";
            }
            std::cout << ss_Y.str()<< std::endl;
        } else {
            // not sure what to do here. Adding 0.1% errors to the perfect predictions?
            for (int i = 0; i < num_samples; ++i) {
                SSE_full += 0.001 * 0.001;
            }
        }
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
        std::cout << "Too few samples (" + std::to_string(num_samples) + ") compared to alleles+covariates (" + std::to_string(num_params_full) + ") in this snarl. Skipping. Note: increasing the minimum number of individuals with -I could help avoiding those issues and get more robust associations in general."<< std::endl;
    } else {
        double numerator = (SSE_reduced - SSE_full) / df_numerator;
        double denominator = SSE_full / df_denominator;
        F_stat = numerator / denominator;
        if (F_stat < 0 && F_stat > -0.00001){
            std::cout << "F statistic is negative but very close to 0 (" + std::to_string(F_stat) + "). Assuming it's 0."<< std::endl;
            F_stat = 0;
        }
        if (F_stat < -0.00001){
            std::cout << "F statistic is negative: " + std::to_string(F_stat) + ". This is concerning, skipping. Recommendation: increase the minimum number of individuals with -I to get more robust associations and avoid issues."<< std::endl;
            std::cout << "SSE_reduced=" + std::to_string(SSE_reduced) +
                        " SSE_full=" + std::to_string(SSE_full) +
                        " df_num=" + std::to_string(df_numerator) +
                        " df_den=" + std::to_string(df_denominator) + " = " + "num_predictors : " + std::to_string(num_predictors) + " num_covariates: " + std::to_string(num_covariates) + " num_samples: " +  std::to_string(num_samples) +
                        " numerator=(" + std::to_string(SSE_reduced) + "-" + std::to_string(SSE_full) + ")/" + std::to_string(df_numerator) +
                        "=" + std::to_string(numerator) +
                        " denominator=" + std::to_string(SSE_full) + "/" + std::to_string(df_denominator) +
                        "=" + std::to_string(denominator) +
                        " F_stat=" + std::to_string(F_stat)<< std::endl;

            std::stringstream ss;
            ss << "X_filtered:\n";
            for (int i = 0; i < X_filtered.rows(); ++i) {
                for (int j = 0; j < X_filtered.cols(); ++j) {
                    ss << X_filtered(i,j) << " ";
                }
                ss << "\n";
            }
            std::cout << ss.str()<< std::endl;

            std::stringstream ss_Y;
            ss_Y << "Y:\n";
            for (int i = 0; i < Y.size(); ++i) {
                ss_Y << Y[i] << " ";
            }
            std::cout << ss_Y.str()<< std::endl;
        }
    }

    // compute a P-value
    boost::math::fisher_f_distribution<double> dist(df_numerator, df_denominator);
    double p_value = 1.0 - boost::math::cdf(dist, F_stat);
    std::cout << "p_value : " << p_value << std::endl;
}

// g++ -std=c++17 -I/usr/include/eigen3 -lboost_math_c99 -o f_stat f_stat.cpp && ./f_stat
