// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Generate random data based on Monte Carlo Simulation
// [[Rcpp::export]]
void graph2dataCpp(arma::mat& D, const arma::mat& A, size_t n=200, size_t iter=15,
                        double val=100.0, double sd=2.0, double prop=0.1, double noise=1.0, bool init=false) {

    // Work with transpose, more efficient considering major column format
    arma::mat A_t = A.t();

    // Get number of columns
    size_t A_ncols = A.n_cols;

    // For each data column
    for (size_t i=0; i < n; i++) {

        // Define container for data column
        arma::vec X(A_ncols);

        if (init) {
            // Take data column of init
            X = D.col(i);
        }
        else {
            // Generate random data column vector
            X = Rcpp::rnorm(A_ncols, val, sd);
        }

        // Run iterations
        for (size_t j=0; j < iter; j++) {
            // Randomize node list
            arma::uvec n_idx = arma::randperm(A_ncols);

            // Iterate through each source node
            for (size_t node=0; node < A_ncols; node++) {
                // Randomize target node list
                arma::uvec targets_tmp = arma::find(A_t.col(n_idx(node)));
                arma::uvec rand_seq = arma::randperm(targets_tmp.n_rows);
                arma::uvec t_idx = targets_tmp(rand_seq);

                // Iterate through target node list
                for (size_t target=0; target < t_idx.n_rows; target++) {
                    // Edge weight, usually 1
                    // MN Has to be std::abs to explicitly apply to floats
                    double P = std::abs( A(n_idx(node), t_idx(target)) );

                    // Updated value for target:
                    // take small fraction of source
                    // add it to large fraction of target
                    // over time values become more and more similar
                    double nval = ( X(n_idx(node)) * prop * P ) +
                        ( X(t_idx(target)) * (1 - prop * P) );

                    // In case of negative associations
                    if ( A(n_idx(node), t_idx(target)) < 0) {
                        double delta = nval - X(t_idx(target));
                        nval = X(t_idx(target)) - delta;
                    }

                    X(t_idx(target)) = nval;
                }
            }
            // Add noise
            X += Rcpp::rnorm(A_ncols, 0, noise);
        }
        // Update resulting data matrix
        D.col(i) = X;
    }
}