#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat mat_prod_prod(arma::mat A, arma::mat B) {
  int I = A.n_cols;
  int K = B.n_cols;
  arma::mat ret(I, K);

  for (int i = 0; i <I; i++) {
    ret.row(i) = arma::prod(1 - (A.col(i) % B.each_col()), 0);
  }
  return ret;
}

// [[Rcpp::export]]
arma::mat A_times_tB(arma::rowvec A, arma::mat B) {
  return (A % B.each_row()).t();
}

// [[Rcpp::export]]
NumericVector get_vec_pos_test(arma::mat P_ik, arma::mat A_ji, arma::mat R_jk, arma::mat Q_jk) {
  arma::mat P_ik_comp = 1 - P_ik;
  arma::mat Q_jk_comp = Q_jk - 1;
  int I = A_ji.n_cols;
  NumericVector vec_pos_test(I);
  for (int i = 0; i < I; i++) {
    arma::mat B = 1 - A_ji.col(i) % R_jk.each_col();
    arma::mat A_and_R = A_times_tB(P_ik.row(i), B);
    arma::mat A_and_Q = A_times_tB(P_ik_comp.row(i), 1 + A_ji.col(i) % Q_jk_comp.each_col());
    vec_pos_test[i] = 1 - arma::prod(arma::prod(A_and_R + A_and_Q));
  }
  return vec_pos_test;
}