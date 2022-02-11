// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Function to perform greedy variable selection based on a correlation matrix
// [[Rcpp::export]]
void lvsCpp(arma::mat &A, arma::mat &RS, const arma::mat &cor_mat, const arma::mat &D, double rs=0.04, int k=5) {

    // Dimensions, cols = variables, rows = samples
    int n = D.n_rows;
    int p = D.n_cols;

    // Expand matrix for getting also the intercept
    colvec C(n, fill::ones);
    arma::mat X = join_horiz(C, D);

    // Main loop
    for (int i=0; i < p; i++) {

        // Indices for set of candidates
        uvec idx_cand = sort_index(arma::abs(cor_mat.col(i)), "descend");

        // Reduced set of indices for set of candidates
        idx_cand = idx_cand(span(0, k - 1));

        // Temporary node list
        vec chosen_tmp(idx_cand.n_rows);

        // Fill with negative values to discriminate
        chosen_tmp.fill(-1);

        // Null model
        colvec coef = solve(X.col(0), D.col(i));
        colvec resid = D.col(i) - X.col(0) * coef;
        double RSS = sum(square(resid));
        // Use version of AIC: -2 * max(logLik) + 2 * p = n * log(RSS / n) + 2 * p
        double AIC = (double) n * log(RSS / (double) n) + 2.0 * X.col(0).n_cols;
        double Rsq = 0.0;

        // Go through all candidates and find the best model
        for (unsigned int j=0; j < idx_cand.n_rows; j++) {

            // Include candidate node
            chosen_tmp(j) = idx_cand(j);

            // Additional set of indices for solver to get also intercept, i.e. adapt index with + 1
            uvec idx_tmp = find(chosen_tmp >= 0);
            uvec idx_intercept = {0};
            uvec idx = join_vert(idx_intercept, idx_cand(idx_tmp) + 1);

            // Linear model X ~ y with intercept
            colvec coef = solve(X.cols(idx), D.col(i));
            colvec resid = D.col(i) - X.cols(idx) * coef;

            // Residual sum of squares: e_sum = sum((y_i - X*coef)^2)
            double RSS = sum(square(resid));

            // Total sum of squares: sum((y_i - mean(y))^2)
            double TSS = sum(square(D.col(i) - (sum(D.col(i)) / n)));

            // Independent parameters are given by the coefficients, including the intercept
            double AIC_tmp = (double) n * log(RSS / n) + 2.0 * ((double) idx.n_rows);
            double Rsq_tmp = 1.0 - (RSS / TSS);

            // Check, whether model gets actually better
            if ( (AIC_tmp < AIC) && (Rsq_tmp - Rsq >= rs) ) {

                // Assign the R^2 of the new node explaining the variation to this degree of the source node to the matrix
                RS(i, idx_cand(j)) = Rsq_tmp - Rsq;

                // Update AIC & Rsq & node list
                AIC = AIC_tmp;
                Rsq = Rsq_tmp;

                // Update resulting adjacency matrix based on sign of correlation
                // Just in case of cor == 0, dont include any edge
                if (cor_mat(i, idx_cand(j)) < 0 ) {
                    A(i, idx_cand(j)) = -1;
                }
                else if (cor_mat(i, idx_cand(j)) > 0) {
                    A(i, idx_cand(j)) = 1;
                }
            }

            // Dump node, if model does not become better
            else {
                chosen_tmp(j) = -1;
            }
        }
    }

    // Convert to undirected graph, but consider negative edges
    // The case that edges will nullify, seems not probable
    A +=  A.t();
    A.replace(2, 1);
    A.replace(-2, -1);
}