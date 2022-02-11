#include <Rcpp.h>

using namespace Rcpp;

// Decision tree based variable selection in C++ / R hybrid
// [[Rcpp::export]]
List tvsRcpp(const DataFrame &D, double rs=0.04, int k=5) {

    // R IMPORTS //
    //  'Import' 'rpart' function
    Environment rpart("package:rpart");                     // Get function from environment
    Function rpart_r = rpart["rpart"];                      // Make it callable from C++

    // 'Import' 'cor' function for method 'spearman' & 'predict' function
    Environment stats("package:stats");
    Function cor_r = stats["cor"];
    Function predict_r = stats["predict"];

    // 'Import' 'paste0_r'
    Environment base("package:base");
    Function paste0_r = base["paste0"];

    // Data frame dimensions
    int p = D.size();

    // Resulting adjacency matrix
    NumericMatrix A = no_init(p, p);
    A.fill(0);

    // R^2 matrix
    NumericMatrix RS = no_init(p, p);
    RS.fill(0);

    // Spearman correlation
    NumericMatrix cor_mat = cor_r(_["x"] = D,
                                  _["use"] = "pairwise.complete.obs",
                                  _["method"] = "spearman");
    // Naming
    rownames(cor_mat) = colnames(D);
    colnames(cor_mat) = colnames(D);
    cor_mat.fill_diag(-2.0);                              // Not include diagonal in further processing

    for (int node=0; node < p; node++) {

        // String representation of nodes
        CharacterVector nodenames = D.names();

        // Finding the node column
        NumericVector targetlist = cor_mat.column(node);

        // Sort Rcpp container, do reverse to get decreasing order
        NumericVector targetlist_sorted = rev(targetlist.sort());
        IntegerVector cand_idx = match(targetlist_sorted, cor_mat.column(node)) - 1;

        // Create positional indices, Note: this starts from 1, so - 1
        IntegerVector subset = seq_len(k) - 1;
        cand_idx = cand_idx[subset];

        // Vector of chosen nodes
        NumericVector chosen(cand_idx.size());

        // R^2 of null model
        double Rsq = 0.0;

        // Get string representation for the variabels to be tested
        String source = nodenames(node);

        // Go through targets
        for (int target=0; target < cand_idx.size(); target++) {

            // Include candidate node
            chosen(target) = cand_idx(target);

            // Already chosen nodes
            // We have to do this, because we update this vector for each iteration
            IntegerVector chosen_idx;

            // Collect already found chosen nodes
            for (int i=0; i < chosen.size(); i++) {
                if (chosen(i) > 0) {
                    chosen_idx.push_back(chosen(i));
                }
            }

            // Get the string representation for the indices
            CharacterVector candidates = nodenames[chosen_idx];

            // If no candidates skip
            if (candidates.size() != 0) {

                String candidates_string = paste0_r(_["..."] = candidates,
                                                    _["collapse"] = "+");

                CharacterVector formula;
                formula.push_back(source);
                formula.push_back("~");
                formula.push_back(candidates_string);

                String formula_str = paste0_r(_["..."] = formula,
                                            _["collapse"] = "");

                List rpart_res = rpart_r(_["formula"] = formula_str,
                                        _["data"] = D);

                NumericVector prediction = predict_r(_["object"] = rpart_res,
                                                    _["newdata"] = D);
                NumericVector D_col = D[node];
                NumericVector Rsq_vec = cor_r(_["x"] = prediction,
                                    _["y"] = D_col,
                                    _["use"] = "pairwise.complete.obs",
                                    _["method"] = "spearman");
                double Rsq_tmp = Rsq_vec[0];
                Rsq_tmp = pow(Rsq_tmp, 2);

                //Check if R^2 gets better
                if (Rsq_tmp - Rsq >= rs) {
                    RS(node, cand_idx(target)) = Rsq_tmp - Rsq;
                    Rsq = Rsq_tmp;

                    if (cor_mat(node, cand_idx(target)) < 0 ) {
                        A(node, cand_idx(target)) = -1;
                    }
                    else if (cor_mat(node, cand_idx(target)) > 0) {
                        A(node, cand_idx(target)) = 1;
                    }
                }

                else {
                    chosen(target) = -1;
                }
            }
        }
    }

    // Conversion to undirected graph should be handeled in R, is faster
    return List::create(Named("adjacency") = A, Named("rs_matrix") = RS);
}