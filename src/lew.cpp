#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

List calculateBetasAndSE(Eigen::Map<Eigen::MatrixXd> XX,
                         Eigen::Map<Eigen::MatrixXd> db_CpG,
                         int trait) {

  Eigen::MatrixXd XtXinv = (XX.transpose() * XX).inverse();
  double XtXinv_se_arg = sqrt(XtXinv(trait, trait));
  int numExplan = XX.cols();

  Eigen::MatrixXd XXproj = XtXinv * XX.transpose();
  Eigen::MatrixXd betas_mat = XXproj * db_CpG;
  Eigen::VectorXd betas = betas_mat.row(trait);

  Eigen::MatrixXd resid_Ys = db_CpG - XX * XXproj * db_CpG;
  Eigen::VectorXd sum_squares_resids = resid_Ys.array().square().colwise().sum();
  Eigen::VectorXd sigmas_square = sum_squares_resids / (db_CpG.rows() - numExplan);
  Eigen::VectorXd se_betas = sigmas_square.array().sqrt() * XtXinv_se_arg;
  Eigen::VectorXd sebetas = ( (1.0 / (XX.rows() - 2)) * sum_squares_resids.array() / ((XX.array() - XX.mean()).array().square().sum() ) ).sqrt();

  return List::create(Named("betas") = betas,
                      Named("se_betas") = se_betas,
                      Named("sebetas") = sebetas);
}

// [[Rcpp::export]]

Eigen::VectorXd solveAndCalculateResiduals(Eigen::MatrixXd& XX, Eigen::MatrixXd& db_CpG) {
  // Calculate coefficients for the full model
  Eigen::MatrixXd b_full = (XX.transpose() * XX).inverse() * XX.transpose() * db_CpG;

  // Calculate sums of squares residuals for the full model
  Eigen::MatrixXd residuals_full = (db_CpG - XX * b_full).array().square().matrix();

  Eigen::VectorXd ss_res_full = residuals_full.colwise().sum();

  return ss_res_full;
}
