// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <queue>

using namespace arma;

// [[Rcpp::export]]
void bfsCpp(arma::mat &S, const arma::mat &G) {

    size_t n = G.n_cols;

    // Init resulting matrix
    // arma::mat S(n, n);
    // S.fill(std::numeric_limits<double>::infinity());
    // S.diag().fill(0);

    // Main loop
    for (size_t i = 0; i < n; i++) {

        // Init queue & mark vector
        std::deque<size_t> Q;
        arma::vec mark = zeros<vec>(n);

        // Set node as discovered
        mark(i) = 1;

        // Include node into queue
        Q.push_back(i);

        while (!Q.empty()) {

            // Get first element of queue
            size_t u = Q.front();
            Q.pop_front();

            // Find all neighbors of current node
            uvec N = find(G.col(u) != 0);

            for (uword ng=0; ng < N.n_rows; ng++) {
                if (S(i, N(ng)) > S(i, u) + 1) {
                    S(i, N(ng)) = S(i, u) + 1;
                }
                if (mark(N(ng)) != 1) {
                    mark(N(ng)) = 1;
                    Q.push_back(N(ng));
                }
            }
        }
    }
}