// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

 //[[Rcpp::export]]
void frgCpp(arma::mat &P, const arma::mat &A, arma::uword iter = 500, double W = 10.0, double L = 10.0,
            double temp_prop = 0.1, double force_prop = 1, double quench_prop = 0.9, double simmering_prop = 0.05) {
    double d, x, y;
    arma::uword n = A.n_cols;
    double area = (W * L);
    double k = force_prop * std::sqrt(area / (double) n);
    double t = W * temp_prop;
    double t_init = t;
    arma::uvec idx = arma::find(A != 0);
    arma::uvec us = arma::floor(idx / n);
    arma::uvec vs = idx - (n  * us);
    arma::mat disp(n, 2);

    // REPULSIVE FORCES
    for (arma::uword i = 0; i < iter; i++) {
        for (arma::uword v = 0; v < n; v++) {
            disp(v, 0) = 0.0;
            disp(v, 1) = 0.0;
            for (arma::uword u = 0; u < n; u++) {
                if (v != u) {
                    // Difference postion between nodes
                    double x = P(v, 0) - P(u, 0);
                    double y = P(v, 1) - P(u, 1);
                    double delta = std::sqrt(x * x + y * y);
                    if (delta < 2.0 * k) {
                        d = k * k / (delta * delta);
                        disp(v, 0) += x * d;
                        disp(v, 1) += y * d;
                    }
                }
            }
        }

        // FORCES OF ATTRACTION
        for (arma::uword e = 0; e < us.n_rows; e++) {
            double x = P(vs(e), 0) - P(us(e), 0);
            double y = P(vs(e), 1) - P(us(e), 1);
            double delta = std::sqrt(x * x + y * y);
            d = (delta * delta) / k * delta;
            disp(vs(e), 0) -= x * d;
            disp(vs(e), 1) -= y * d;
            disp(us(e), 0) += x * d;
            disp(us(e), 1) += y * d;
        }

        // LIMITING
        for (arma::uword v = 0; v < n; v++) {
            double dx = disp(v, 0);
            double dy = disp(v, 1);
            double disp = std::sqrt(dx * dx + dy * dy);
            if (disp != 0.0) {
                // Limit displacement to be at maximum t
                d = std::min(disp, t) / disp;
                x = P(v, 0) + dx * d;
                y = P(v, 1) + dy * d;

                // Avoid points to be outside of the frame
                double px = (W / 2) - (y - L / 2) + (W / 2);
                double py = (L / 2) - (x - W / 2) + (L / 2);
                P(v, 0) = std::min(W / 2.0 - py + L / 2 + W / 2, std::max(-1 * (W / 2.0 - py + L / 2 + W / 2), x));
                P(v, 1) = std::min(L / 2.0 - px + W / 2 + L / 2, std::max(-1 * (L / 2.0 - px + W / 2 + L / 2), y));
            }
        }

        // COOLING
        if (t >= t_init * simmering_prop) {
            t *= pow(quench_prop, i);
        }
    }
}