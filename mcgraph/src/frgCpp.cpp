// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

 //[[Rcpp::export]]
void frgCpp(arma::mat &P, const arma::mat &A, arma::uword iter = 500, double W = 10.0, double L = 10.0,
            double temp_prop = 0.1, double force_prop = 0.5, double quench_prop = 0.9, double simmering_prop = 0.05) {
    // shrink initial values
    P * 1e-6;
    double area, k, t, t_init, d, delta_x, delta_y, fr, fa, d_disp, width_side, length_side;
    arma::uword n = A.n_cols;
    area = (W * L);
    // k is the radius
    k = force_prop * std::sqrt(area / (double) n);
    t = W * temp_prop;
    t_init = t;
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
                    // difference postion between nodes
                    delta_x = P(v, 0) - P(u, 0);
                    delta_y = P(v, 1) - P(u, 1);
                    // Euclidean distance
                    d = std::sqrt(delta_x * delta_x + delta_y * delta_y);
                    // repulse only, if nodes are close to each other
                    if (2.0 * k >= d) {
                        //double fr = 1;
                        fr = (k * k / d);
                        disp(v, 0) += (delta_x / d) * fr;
                        disp(v, 1) += (delta_y / d) * fr;
                    }
                }
            }
        }

        // FORCES OF ATTRACTION
        for (arma::uword e = 0; e < us.n_rows; e++) {
            delta_x = P(vs(e), 0) - P(us(e), 0);
            delta_y = P(vs(e), 1) - P(us(e), 1);
            d = std::sqrt(delta_x * delta_x + delta_y * delta_y);
            //double fa = 1;
            fa = (d * d) / k;
            disp(vs(e), 0) -= (delta_x / d) * fa;
            disp(vs(e), 1) -= (delta_y / d) * fa;
            disp(us(e), 0) += (delta_x / d) * fa;
            disp(us(e), 1) += (delta_y / d) * fa;
        }

        // LIMITING
        for (arma::uword v = 0; v < n; v++) {
            // displacement distance
            d_disp = std::sqrt(disp(v, 0) * disp(v, 0) + disp(v, 1) * disp(v, 1));
            // just limit for actual displacement
            if (d_disp != 0.0) {
                // limit displacement to be at maximum t
                P(v, 0) += (disp(v, 0) / d_disp) * std::min(d_disp, t);
                P(v, 1) += (disp(v, 1) / d_disp) * std::min(d_disp, t);
                // avoid points to be outside of the frame, model the frame as a circle
                width_side = (W / 2) - (P(v, 1) - L / 2) + (W / 2);
                length_side = (L / 2) - (P(v, 0) - W / 2) + (L / 2);
                P(v, 0) = std::min(W / 2.0 - length_side + L / 2 + W / 2, std::max(-(W / 2.0 - length_side + L / 2 + W / 2), P(v, 0)));
                P(v, 1) = std::min(L / 2.0 - width_side + W / 2 + L / 2, std::max(-(L / 2.0 - width_side + W / 2 + L / 2), P(v, 1)));
            }
        }

        // COOLING
        if (t >= t_init * simmering_prop) {
            t *= pow(quench_prop, i);
        }
    }
}