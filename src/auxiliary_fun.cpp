// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;
const double log2pi = std::log(2.0 * M_PI);

//' Multivariate normal density function
//'
//' This function returns the log-density for a multivariate Gaussian distribution.
//' The data must be imputed as a matrix, using e.g., \code{as.matrix}, with each row
//' representing an observation.
//'
//' @param x matrix of observations
//' @param mean mean vector
//' @param sigma positive definite covariance matrix
//' @param logd logical; whether log-density should be returned (default to \code{FALSE})
//' @return density or log-density of the \code{nrow(x)} sample
//' @keywords internal
//' @export
// [[Rcpp::export(.dmvnorm_arma)]]
arma::vec dmvnorm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma,
                      bool logd = false) {
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
  }

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

//' Distance matrix with geometric anisotropy
//'
//' The function computes the distance between locations, with geometric anisotropy.
//' The parametrization assumes there is a scale parameter, so that \code{scale}
//' is the distortion for the second component only. The angle \code{rho} must lie in
//' \eqn{[-\pi/2, \pi/2]}.
//'
//' @param loc a \code{d} by 2 matrix of locations giving the coordinates of a site per row.
//' @param scale numeric vector of length 1, greater than 1.
//' @param rho angle for the anisotropy, must be larger than \eqn{\pi/2} in modulus.
//' @return a \code{d} by \code{d} square matrix of pairwise distance
//' @export
// [[Rcpp::export]]
arma::mat distg(arma::mat loc, NumericVector scale, NumericVector rho){
  int d = loc.n_rows;
  int m = loc.n_cols;
  arma::mat aniso(2, 2);
  if(loc.n_cols > 2){
    stop("Invalid location matrix; only geometric anisotropy for bivariate is supported");
  }
  if(rho.size() > 1 || scale.size() > 1){
    stop("Invalid length for `scale` or `rho`");
  }
  if(std::abs(rho[0]) > M_PI){
    stop("Invalid `rho` argument: angle must be in [-pi/2, pi/2]");
  }
  if(scale[0] < 1){
    stop("Scale parameter should be larger than 1 for identifiability");
  }
    aniso(0,0) = cos(rho)[0];
      aniso(0,1) = sin(rho)[0];
      aniso(1,0) = -scale[0]*aniso(0,1);
      aniso(1,1) = scale[0]*aniso(0,0);
      arma::mat distmat(d, d);
      distmat.zeros();
      arma::vec dist_ij(loc.n_cols);
      for(int i=0; i<d-1; i++){
        for(int j=i+1; j < d; j++){
          for(int k=0; k < m; k++){
            dist_ij(k) = loc(i,k) - loc(j,k);
          }
          distmat(i,j) = arma::norm(aniso * dist_ij);
          distmat(j,i) = distmat(i,j);
        }
      }
      return distmat;
}



