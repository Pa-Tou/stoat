// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::plugins(openmp)]]

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT

#include <stdio.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <set>
#include <vector>
#include <iostream>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppGSL.h>

#include "Eigen/Cholesky"
#include "Eigen/Core"
#include "gsl/gsl_cdf.h"
#include <gsl/gsl_cdf.h>

#define DECLARE_EIGEN_VECTOR(v, v_e) Eigen::Map<Eigen::VectorXd> v_e((v).data.data(), (v).data.size())
#define DECLARE_EIGEN_CONST_VECTOR(v, v_e) Eigen::Map<const Eigen::VectorXd> v_e((v).data.data(), (v).data.size())
#define DECLARE_EIGEN_MATRIX(m, m_e) Eigen::Map<Eigen::MatrixXd> m_e((m).data.data(), (m).rows, (m).cols)
#define DECLARE_EIGEN_CONST_MATRIX(m, m_e) Eigen::Map<const Eigen::MatrixXd> m_e((m).data.data(), (m).rows, (m).cols)

// ======================= Vector Class =========================
class Vector {
  public:
  std::vector<double> data;

  Vector() {}
  Vector(int n);
  Vector(int n, double val);

  double& operator[](int i);
  double operator[](int i) const;
  int Length() const;
  void Dimension(int n);
  void Dimension(int n, double val);
  void Fill(double val);
  double Sum() const;
  double Average() const;
  double Min() const;
  double Max() const;
};

Vector::Vector(int n) { Dimension(n); }
Vector::Vector(int n, double val) { Dimension(n, val); }
void Vector::Dimension(int n) { data.resize(n); }
void Vector::Dimension(int n, double val) {
  data.resize(n);
  Fill(val);
}
void Vector::Fill(double val) { std::fill(data.begin(), data.end(), val); }
double Vector::Sum() const {
  return std::accumulate(data.begin(), data.end(), 0.0);
}
double Vector::Average() const {
  return data.empty() ? 0.0 : Sum() / data.size();
}
double Vector::Min() const {
  return *std::min_element(data.begin(), data.end());
}
double Vector::Max() const {
  return *std::max_element(data.begin(), data.end());
}

int Vector::Length() const { return data.size(); }
double& Vector::operator[](int i) { return data[i]; }
double Vector::operator[](int i) const { return data[i]; }

// ======================= Matrix Class =========================
class Matrix {
  public:
  int rows, cols;
  std::vector<double> data;
  std::vector<std::string> colLabel;

  Matrix();
  explicit Matrix(int nr, int nc);
  Matrix(const Matrix& m);
  Matrix& operator=(const Matrix& m);

  double& operator()(int r, int c) { return data[r + c * rows]; }
  double operator()(int r, int c) const { return data[r + c * rows]; }

  void Dimension(int nr, int nc);
  void Dimension(int nr, int nc, double val);
  void DimensionQuick(int nr, int nc);
  void Reserve(int nr, int nc);
  void Fill(double val);
  double Min() const;
  double Max() const;
  int Length() const;

  void Product(const Matrix& in1, const Matrix& in2);
  void Transpose(const Matrix& old);
  Matrix& Multiply(double s);
  int RemoveByRowIndex(const std::vector<int>& rowIndexToRemove);
  Matrix& StackRight(const Matrix& m);
};

Matrix::Matrix() {
  rows = 0;
  cols = 0;
}

Matrix::Matrix(int nr, int nc) : rows(nr), cols(nc), data(nr * nc) {}

Matrix::Matrix(const Matrix& m)
    : rows(m.rows), cols(m.cols), data(m.data), colLabel(m.colLabel) {}

Matrix& Matrix::operator=(const Matrix& m) {
  rows = m.rows;
  cols = m.cols;
  data = m.data;
  colLabel = m.colLabel;
  return *this;
}

void Matrix::Dimension(int nr, int nc) {
  if (nr == rows && nc == cols) return;
  std::vector<double> newData(nr * nc);
  for (int i = 0; i < nr && i < rows; ++i)
    for (int j = 0; j < nc && j < cols; ++j)
      newData[i + j * nr] = data[i + j * rows];
  rows = nr;
  cols = nc;
  std::swap(data, newData);
  colLabel.resize(nc);
}

int Matrix::Length() const { return data.size(); }

void Matrix::Dimension(int nr, int nc, double val) {
  DimensionQuick(nr, nc);
  Fill(val);
}

void Matrix::DimensionQuick(int nr, int nc) {
  rows = nr;
  cols = nc;
  data.resize(nr * nc);
  colLabel.resize(nc);
}

void Matrix::Reserve(int nr, int nc) {
  data.reserve(nr * nc);
  colLabel.reserve(nc);
}

void Matrix::Fill(double val) {
  std::fill(data.begin(), data.end(), val);
}

double Matrix::Min() const {
  return *std::min_element(data.begin(), data.end());
}

double Matrix::Max() const {
  return *std::max_element(data.begin(), data.end());
}

void Matrix::Product(const Matrix& in1, const Matrix& in2) {
  DECLARE_EIGEN_CONST_MATRIX(in1, e1);
  DECLARE_EIGEN_CONST_MATRIX(in2, e2);
  DimensionQuick(in1.rows, in2.cols);
  DECLARE_EIGEN_MATRIX((*this), out);
  out = e1 * e2;
}

void Matrix::Transpose(const Matrix& old) {
  data.resize(old.data.size());
  rows = old.cols;
  cols = old.rows;
  DECLARE_EIGEN_CONST_MATRIX(old, eOld);
  DECLARE_EIGEN_MATRIX((*this), eNew);
  eNew = eOld.transpose();
}

Matrix& Matrix::Multiply(double s) {
  for (double& val : data) val *= s;
  return *this;
}

int Matrix::RemoveByRowIndex(const std::vector<int>& rowIndexToRemove) {
  int idx = 0;
  std::set<int> idxSet(rowIndexToRemove.begin(), rowIndexToRemove.end());
  for (int j = 0; j < cols; ++j) {
    for (int i = 0; i < rows; ++i) {
      if (idxSet.count(i)) {
        continue;
      }
      data[idx++] = (*this)(i, j);
    }
  }
  rows -= idxSet.size();
  data.resize(rows * cols);
  return idxSet.size();
}

Matrix& Matrix::StackRight(const Matrix& m) {
  data.insert(data.end(), m.data.begin(), m.data.end());
  cols += m.cols;
  colLabel.insert(colLabel.end(), m.colLabel.begin(), m.colLabel.end());
  return *this;
}

// ==================== LinearRegression Class ===================
struct LinearRegression {

  Matrix XtXinv, covB;
  Vector predict, residuals, B, pValue;
  double sigma2;

  LinearRegression() : sigma2(0.){};
  bool FitLinearModel(const Matrix& X, const Vector& y);
  Vector& GetAsyPvalue();
  Vector& GetCovEst() { return this->B; };  // (X'X)^{-1} X'Y
  Matrix& GetCovB() { return this->covB; };
  Vector& GetPredicted() { return this->predict; };
  Vector& GetResiduals() { return this->residuals; };
  double GetSigma2() const { return this->sigma2; };
  bool calculateHatMatrix(Matrix& X, Matrix* out);
  bool calculateResidualMatrix(Matrix& X, Matrix* out);
};

bool LinearRegression::FitLinearModel(const Matrix& X, const Vector& y) {
  XtXinv.Dimension(X.cols, X.cols);
  DECLARE_EIGEN_CONST_MATRIX(X, X_e);
  DECLARE_EIGEN_MATRIX(XtXinv, XtXinv_e);
  XtXinv_e = (X_e.transpose() * X_e).llt().solve(Eigen::MatrixXd::Identity(X.cols, X.cols));

  B.Dimension(X.cols, 1);
  DECLARE_EIGEN_VECTOR(B, B_e);
  DECLARE_EIGEN_CONST_VECTOR(y, y_e);
  B_e = XtXinv_e * X_e.transpose() * y_e;

  predict.Dimension(X.rows, 1);
  residuals.Dimension(X.rows, 1);
  DECLARE_EIGEN_VECTOR(predict, pred_e);
  DECLARE_EIGEN_VECTOR(residuals, resid_e);
  pred_e = X_e * B_e;
  resid_e = y_e - pred_e;

  sigma2 = resid_e.squaredNorm() / y.Length();

  covB.Dimension(X.cols, X.cols);
  DECLARE_EIGEN_MATRIX(covB, covB_e);
  covB_e = XtXinv_e * sigma2;

  return true;
}

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::export]]
double my_gsl_chisq_Q(double x, double df) {
  return gsl_cdf_chisq_Q(x, df);
}

Vector& LinearRegression::GetAsyPvalue() {
  int numCov = B.Length();
  pValue.Dimension(numCov);
  for (int i = 0; i < numCov; ++i) {
    double se = sqrt(covB(i, i));
    double t_stat = B[i] / se;

    double Zstat = t_stat * t_stat;
    pValue[i] = my_gsl_chisq_Q(Zstat, 1.0); // why degree of freedom = 1.0 ?
  }

  return pValue;
}

bool LinearRegression::calculateHatMatrix(Matrix& X, Matrix* out) {
  DECLARE_EIGEN_CONST_MATRIX(X, X_e);
  XtXinv.Dimension(X.cols, X.cols);
  DECLARE_EIGEN_MATRIX(XtXinv, XtXinv_e);
  XtXinv_e = (X_e.transpose() * X_e).llt().solve(Eigen::MatrixXd::Identity(X.cols, X.cols));
  DECLARE_EIGEN_MATRIX((*out), out_e);
  out_e = X_e * XtXinv_e * X_e.transpose();
  return true;
}

bool LinearRegression::calculateResidualMatrix(Matrix& X, Matrix* out) {
  if (!calculateHatMatrix(X, out)) return false;
  for (int i = 0; i < out->rows; ++i)
    for (int j = 0; j < out->cols; ++j)
      (*out)(i, j) = (i == j) ? 1.0 - (*out)(i, j) : -(*out)(i, j);
  return true;
}

// [[Rcpp::export]]
Rcpp::List cpp_linear_regression_rvtest(Rcpp::NumericMatrix Xr, Rcpp::NumericVector yr) {
    int numSamples   = Xr.nrow();
    int numVariables = Xr.ncol();

    // Convert Rcpp matrix to your C++ Matrix
    ::Matrix X(numSamples, numVariables); // explicitly global
    for (int i = 0; i < numSamples; i++) {
        for (int j = 0; j < numVariables; j++) {
            X(i, j) = Xr(i, j);
        }
    }

    // Convert phenotype vector
    ::Vector y; // explicitly global
    y.Dimension(numSamples);
    for (int i = 0; i < numSamples; i++) {
        y[i] = yr[i];
    }

    LinearRegression lr;
    bool success = lr.FitLinearModel(X, y);

    if (!success) {
        Rcpp::stop("Linear regression failed.");
    }

    ::Vector &coef = lr.GetCovEst();
    ::Vector &pval = lr.GetAsyPvalue();

    Rcpp::NumericVector coef_out(numVariables);
    Rcpp::NumericVector pval_out(numVariables);
    for (int i = 0; i < numVariables; i++) {
        coef_out[i] = coef[i];
        pval_out[i] = pval[i];
    }

    return Rcpp::List::create(
        Rcpp::_["coefficients"] = coef_out,
        Rcpp::_["p_values"] = pval_out
    );
}

// ======================= Test Main ===========================
int main() {
  const int numSamples = 9;
  const int numVariables = 3;  // Intercept + 2 SNP variables
  
  Matrix X(numSamples, numVariables);
  
  // Row-wise assignment: [intercept, df[i][0], df[i][1]]
  X(0, 0) = 1.0; X(0, 1) = 1.0; X(0, 2) = 0.0;
  X(1, 0) = 1.0; X(1, 1) = 1.0; X(1, 2) = 0.0;
  X(2, 0) = 1.0; X(2, 1) = 1.0; X(2, 2) = 0.0;
  X(3, 0) = 1.0; X(3, 1) = 1.0; X(3, 2) = 0.0;
  X(4, 0) = 1.0; X(4, 1) = 1.0; X(4, 2) = 0.0;
  X(5, 0) = 1.0; X(5, 1) = 1.0; X(5, 2) = 0.0;
  X(6, 0) = 1.0; X(6, 1) = 1.0; X(6, 2) = 0.0;
  X(7, 0) = 1.0; X(7, 1) = 0.0; X(7, 2) = 1.0;
  X(8, 0) = 1.0; X(8, 1) = 0.0; X(8, 2) = 0.0;
  
  // Phenotype vector
  Vector y;
  y.Dimension(numSamples);
  y[0] = 4.5;
  y[1] = 7.0;
  y[2] = 9.2;
  y[3] = 10.9;
  y[4] = 13.0;
  y[5] = 14.0;
  y[6] = 11.0;
  y[7] = 15.0;
  y[8] = 16.0;  

  LinearRegression lr;
  if (lr.FitLinearModel(X, y)) {
    Vector& coef = lr.GetCovEst();
    printf("Coefficients:\n");
    for (int i = 0; i < coef.Length(); i++) {
      printf("  B[%d] = %f\n", i, coef[i]);
    }

    Vector& pval = lr.GetAsyPvalue();
    printf("P-values:\n");
    for (int i = 0; i < pval.Length(); i++) {
      printf("  p[%d] = %.5g\n", i, pval[i]);
    }
  } else {
    printf("Linear regression failed.\n");
  }

  return 0;
}

// g++ linear_regression_rvtest.cpp -std=c++17 -I/usr/local/include/eigen3 -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm -o linear_regression_rvtest
// ./linear_regression_rvtest
