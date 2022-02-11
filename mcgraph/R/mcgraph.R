#' @title Create random graphs and data for those graphs
#'
#' @description `mcg.graph` creates random directed graphs and simulated data for them from undirected graphs, or just data for directed graphs.
#'
#' @details This function allows you to create random graphs, either directed or undirected and data belonging
#' to the actual graph structure based on Monte Carlo simulations. If the graph is already directed just the data are generated for this graph.
#' @param U adjacency matrix for an undirected graph.
#' @param n number of observations to generate for the created directed graph type type of graph.
#' @param input number of random input nodes, if node labels are given they are taken as input nodes.
#' @return An mcgraph graph object with list item data (simulated data for the directed graph).
#' @author Detlef Groth <email: dgroth@uni-potsdam.de>
#' @keywords Monte Carlo sampling
#' @examples
#' band=mcg.band(16)
#' mcb=mcg.graph(band,n=200,input='A1')
#' plot(mcb,layout='circle')
#' @export

mcg.graph <- function (U,n=200,input=2) {
    if (!is.null(attr(U,"mode")) && attr(U,"mode") == "undirected") {
        D=mcg.u2d(U,input=input)
    } else {
        D=U
    }
    G=D
    attr(G,"data")=mcg.graph2data(G,n=n)
    return(G)
}

# TODO
# * create graph on your own
# * layout as attribute
# * is.mcgraph
# * as.matrix remove attributes
# * corrplot with colors='bluered' 'grey'
# * corrplot with pch.lower = 0 -> numbers
# * negative correleations
# * fix for layout mds two nodes on same point (add noise automatically)
# * graph2data with undirected graph
# * average R with multiple round correlations?
# * negative R's 200-node.val
# * as.matrix -> remove attributes
# * vignette
# * mcg.corrplot
# https://www.r-bloggers.com/sinew-a-r-package-to-create-self-populating-roxygen2-skeletons/

#' @title Create a new `mcgraph` object
#' @description Based on a given adjacency matrix a new mcgraph object is created.
#' @param A Input adjacency matrix.
#' @param type Custom type name for the graph. Default: 'custom'
#' @return A mcgraph object for the given adjacency matrix.
#' @details The function mcg.new generates for the given adjacency matrix a mcgraph object.
#'      If the adjacency matrix is has symmetric upper and lower triangles, the graph is undirected,
#'      otherwise the function will return a directed graph.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  G=matrix(0,nrow=6,ncol=6)
#'  rownames(G)=colnames(G)=LETTERS[1:6]
#'  G['A','C']=1
#'  G['B','C']=1
#'  G['C','D']=1
#'  G['D','E']=1
#'  G['D','F']=1
#'  G['E','F']=1
#'  G=mcg.new(G,type="test")
#'  plot(G)
#'  }
#' }
#' @rdname mcg.new
#' @export

mcg.new <- function (A, type="custom") {
    # Convert in case of e.g. sparse matrix
    if (!is.matrix(A)) {
        A <- as.matrix(A)
    }
    # MN Include squared dims check
    if (nrow(A) != ncol(A)) {
        stop("Matrix is not squared.")
    }
    # MN The new check is 5x faster and takes less memory
    #if (identical(A[lower.tri(A)], t(A)[lower.tri(A)])) {
    if (all(A == t(A))) {
        attr(A,"mode") <- "undirected"
    } else {
       attr(A,"mode") <- "directed"
    }
    # If no dimnames, include them
    if (is.null(colnames(A)) || is.null(rownames(A))) {
        colnames(A) <- rownames(A) <- mcg.autonames(LETTERS, ncol(A))
    }
    attr(A,"type") <- type
    class(A) <- "mcgraph"
    return(A)
}

#' @title Using hard correlation thresholding for graph construction
#' @description mcg.ct - creates a new mcgraph object from a given correlation matrix or using raw data by a simple R-square thresholding mechanism
#' @param x, either a correlation matrix or a data frame / matrix with raw data,
#' @param type custom type name for the graph, Default: 'ct'
#' @param rs threshold for the R-square value. Default: 0.04
#' @param is.squared is the matrix already squared, if not the function will scale the matrix before applying the threshold.
#' @param method for raw data which method for measuring correlation stringth should be used, default: 'pearson'
#' @return an mcgraph object for the given data or correlation  matrix.
#' @details Nodes in the resulting network will be connected, if their pairwise correlations will exceed the given threshold. Negative associations
#'          are set by the negative sign of the respective correlation value.
#' @examples
#'  G=matrix(0,nrow=6,ncol=6)
#'  rownames(G)=colnames(G)=LETTERS[1:6]
#'  G['A','C']=1
#'  G['B','C']=1
#'  G['C','D']=1
#'  G['D','E']=1
#'  G['D','F']=1
#'  G['E','F']=1
#'  G=mcg.new(G,type="test")
#'  G.data=mcg.graph2data(G,n=200)
#'  G2=mcg.ct(cor(t(G.data)))
#'  par(mfrow=c(1,2))
#'  plot(G,main="Real Graph",layout='sam')
#'  plot(G2,main="Predicted Graph",layout='sam')
#'  data(swiss)
#'  plot(mcg.ct(swiss,method="spearman",rs=0.1),layout="star")
#' @rdname mcg.ct
#' @export

mcg.ct <- function (x,type="ct",rs=0.04,is.squared=FALSE,method="pearson") {
    # MN: In case of p == n and if rownames of data is NULL this will not work.
    # if (nrow(x) != ncol(x) || !all(rownames(x)==colnames(x))) {
    #     C=cor(x,method=method,use="pairwise.complete.obs")
    #     if (method == "kendall") {
    #         is.squared=TRUE
    #     } else {
    #         is.squared=FALSE
    #     }
    # } else {
    #     C=x
    # }
    # if (nrow(C) != ncol(C) || !identical(C[lower.tri(C)], t(C)[lower.tri(C)])) {
    #     stop("x must be a symmetric correlation matrix")
    # }
    # MN Check whether the matrix is squared and the first of its values are only drawn from the set {-1, 0, 1}, i.e. whether its an adjacency matrix
    # Weighted graphs are not considered
    # if (ncol(x) == nrow(x) && all(unique(x[1:2,]) %in% c(0, 1, -1))) {
    #     stop("x must be a correlation matrix or a data frame/matrix with raw data, not an adjacency matrix")
    # }
    # MN: If the diag of cor matrix is editted, this will get unnoticed
    # TODO: Find better solution
    if (!all(diag(as.matrix(x)) == 1) & nrow(x) == ncol(x)) {
        stop("x must be a correlation matrix or a data frame/matrix with raw data, not an adjacency matrix")
    }
    # MN Hence, if the matrix is symmetric and square, it must be a correlation matrix
    #if (identical(x[upper.tri(x)], x[lower.tri(x)]) && ncol(x) == nrow(x)) {
    if (all(as.vector(x) == as.vector(t(x))) && ncol(x) == nrow(x)) {
        C <- x
    }
    # MN x contains raw data
    else {
        C <- cor(x, method=method, use="pairwise.complete.obs")
        if (method == "kendall") {
            is.squared=TRUE
        } else {
            is.squared=FALSE
        }
    }
    if (!is.squared) {
        S=C^2
    } else {
        S=abs(C)
    }
    A=S
    A[A > rs]=1
    A[A <=rs]=0
    A[A==1 & C < 0]= -1
    diag(A)=0
    # MN Variables must have colnames
    if (is.null(colnames(x))) {
        colnames(x) <- mcg.autonames(LETTERS, ncol(x))
    }
    rownames(A) <- colnames(A) <- rownames(C) <- colnames(C) <- colnames(x)
    attr(A,"mode")="undirected"
    attr(A,"type")=type
    attr(A,"r.squared")=C
    class(A)="mcgraph"
    return(A)
}

#' @title Forward variable selection for graph construction
#' @description For each variable, first the `k` highest correlated variables are preselected. Variables are included
#'              in the model, if (1) Aikes Information Crition (AIC) becomes smaller and (2) the difference of
#'              Pearson's coefficient of determination (R-square value) becomes larger compared to not including the variable.
#'
#'              \code{NA} entries are internally dealt with by just removing the according rows. Hence, possible imputation should
#'              be performed beforehand.
#' @param d Data frame or matrix with rows including data values and variables in columns. Using only numbers as column names is not allowed.
#' @param rs Threshold for the R-square value. Default: 0.02
#' @param k Maximal number of correlated variables to consider. Default: 5
#' @param output Returned data type, either a mcgraph object, a matrix (adjacency matrix) or a list (adjacency list). Default: 'mcgraph'
#' @param code Should the R or the C++ version be used? Default: 'C++'
#' @param method Method used to calculate the values of the correlation matrix. Default: 'pearson'
#' @return A `mcgraph` matrix object, alternatively a matrix or a list depending on the value of the argument 'output'.
#' @examples
#' data(swiss)
#'
#' Gswiss <- mcg.lvs(swiss, rs=0.1, output='mcgraph')
#' plot(Gswiss, main="Predicted Graph", layout='frg')
#' plot(Gswiss, main="Predicted Graph", layout='star')
#'
#' ang <- mcg.angie(nodes=26, edges=40)
#' data <- mcg.graph2data(ang)
#' plot(ang, layout="frg")
#' @rdname mcg.lvs
#' @author: Masiar Novine <email: masiar.novine@gmail.com>
#' @export

# MN Changes: I removed the cor.squared argument, and work with the absolute correlation values instead
mcg.lvs <- function(d, rs=0.04, k=5, output="mcgraph", code="C++", method="pearson") {
    if (k > ncol(d)) {
        k <- ncol(d)
        warning("Subset of neighbor nodes cannot be equal or larger than number of nodes of graph")
    }
    # MN Variables must have colnames
    if (is.null(colnames(d))) {
        colnames(d) <- mcg.autonames(LETTERS, ncol(d))
    }
    # MN Handle NAs, very strict, just removes the row(s) containing NA
    if (any(is.na(d))) {
        d <- na.omit(d)
    }
    if (code == "C++") {
        # MN 'd' must be S3 matrix
        if (!is.matrix(d)) {
            d <- as.matrix(d)
        }
        # MN Call by reference, init results beforehand, avoid copies
        adj <- matrix(0, nrow=ncol(d), ncol=ncol(d))
        rsqs <- matrix(0, nrow=ncol(d), ncol=ncol(d))
        # TODO: (1) Implement Spearman correlation in C++, 'cor' in R has large overhead,
        # is slower and takes much more memory
        cor.mt <- cor(d, method=method, use="pairwise.complete.obs")
        diag(cor.mt) <- 0
        lvsCpp(adj, rsqs, cor.mt, d, rs, k)
        rownames(adj) <- colnames(adj) <- colnames(d)
        rownames(rsqs) <- colnames(rsqs) <- colnames(d)
    } else if (code == "R") {
        if (!is.data.frame(d)) {
            if (is.matrix(d)) {
                d <- as.data.frame(d)
            } else {
                stop("Error: Input d must be a data frame or matrix!")
            }
        }
        # MN Omit 'cor.squared' argument completely
        cor_d <- cor(d, method=method, use="pairwise.complete.obs")
        # MN Exclude diagonal elements from further calculation
        diag(cor_d) <- 0
        chosen <- list()
        if (any(grepl("^[0-9]", colnames(d)))) {
            stop("Using numbers as column names is not allowed.")
        }
        if (output == "mcgraph" || output == "matrix") {
            adj <- matrix(0, nrow=ncol(d), ncol=ncol(d))
            rownames(adj) <- colnames(adj) <- colnames(d)
            # DG: adding rsquare values
            rsqs <- adj
        }
        for (node in colnames(d)) {
            chosen[node] <- c()
            # MN Use absolute correlation values
            cand <- names(sort(rank(abs(cor_d[node, ])), decreasing=TRUE))[1:(k + 1)]
            # MN Can be omitted, because of zero diagonal
            # cand <- cand[-which(cand == node)]
            model <- lm(formula(paste(c(node, "~", "0"), collapse="")), d)
            aic <- AIC(model)
            # MN Changed, concentrate in one line
            rsq <- summary.lm(model)$adj.r.squared
            for (target in cand) {
                model_tmp <- lm(formula(paste(c(node, "~", paste(c(chosen[[node]], target), collapse="+")), collapse="")), d)
                aic_tmp <- AIC(model_tmp)
                # MN Changed, concentrate in one line
                rsq_tmp <- summary.lm(model_tmp)$adj.r.squared
                if ((aic_tmp < aic && rsq_tmp - rsq >= rs)) {
                    chosen[[node]] <- append(chosen[[node]], target)
                    aic <- aic_tmp
                    if (output == "mcgraph" || output == "matrix") {
                        rsqs[node, target] <- rsq_tmp - rsq
                        if (cor_d[node, target] > 0) {
                            adj[node, target] <- 1
                        }
                        else if (cor_d[node, target] < 0) {
                            adj[node, target] <- -1
                        }
                    }
                    rsq <- rsq_tmp
                }
            }
        }
    } else {
        stop("Error: 'code' must be either 'R' or 'C++'!")
    }
    if (output == "mcgraph") {
        # Already converted in C++ code
        if (code == "R") {
            adj <- adj + t(adj)
            adj[which(adj == 2)] <- 1
            adj[which(adj == -2)] <- -1
        }
        G=mcg.new(adj)
        attr(G,"r.squared") <- rsqs
        return(G)
    }
    else if (output == "matrix") {
        if (code == "R") {
            adj <- adj + t(adj)
            adj[which(adj == 2)] <- 1
            adj[which(adj == -2)] <- -1
        }
        return(adj)
    }
    else if (output == "list") {
        if (code == "R") {
            return(chosen)
        } else {
            # TODO: loop over adj and convert it to a list
            stop("Error: Output 'list' is currently only available for 'code' 'R'")
        }
    }
}

#' @title Variable selection bases on regression trees for graph construction
#' @description The function mcg.rpart - creates a new mcgraph object from
#'        given data using regression trees. The algorithm will evaluate
#'        for every variable the ten variables with the highest absolute
#'        Spearman correlation for their predictive power. Variables which
#'        increase the predictive power, emasured using R-square by more than the threshold, default 0.2, will be connected with edges.
#' @param x data frame or matrix with numerical daattr(Gswiss,"r.squared"),3ta.
#' @param type custom type name for the graph, Default: 'rpart'
#' @param rs threshold for the R-square value. Default: 0.04
#' @param k maximal number of correlated variables to consider, default: 10
#' @param keep what data should be returned, either mcgraph or list
#'        (adjacency list) or matrix (adjacency matrix), list might be useful for
#'        for many columns, default: 'mcgraph'
#' @return an mcgraph object for the data, or a adjacency matrix or list depending on the keep argument,
#'      in the djacency matrix a value of minus indicate negative associations
#'      between the nodes, a value of +1 positive associations.
#' @import rpart
#' @examples
#'  data(swiss)
#'  Gswiss=mcg.rpart(swiss,rs=0.1)
#'  plot(Gswiss,main="Predicted Graph",layout='sam')
#' @rdname mcg.rpart
#' @author Detlef Groth, University of Potsdam
#' @export

mcg.rpart <- function (x,type="rpart",rs=0.04,k=10,keep="mcgraph") {
    rs.threshold=rs
    # TODO: Mutual information as alternative
    C=cor(x,method="spearman",use="pairwise.complete.obs")
    A=matrix(0,nrow=ncol(x),ncol=ncol(x))
    rownames(A)=colnames(A)=colnames(x)
    if (k > ncol(x)) {
        k=ncol(x)
    }
    if (is.matrix(x)) {
        x=as.data.frame(x)
    }
    for (i in 1:ncol(x)) {
        idx=order(abs(C[i,]),decreasing=TRUE)[2:k]
         rp=rpart::rpart(formula(paste(colnames(x)[i], "~",
                                paste( colnames(x)[idx],collapse="+"))),
                                data=x)
        # check all important variables
        vimps=names(rp$variable.importance)
        rp1=rpart::rpart(formula(paste(colnames(x)[i], "~",vimps[1])),data=x)
        rs=cor(predict(rp1,newdata=x),x[,i],use="pairwise.complete.obs",method="spearman")^2
        #return()
        if (rs < rs.threshold) {
            break
        }
        A[i,vimps[1]]=rs

        rss=c(rs)
        #next
        for (vi in 2:length(vimps)) {
            rpi=rpart::rpart(formula(paste(colnames(x)[i], "~",
                                paste( vimps[1:vi],collapse="+"))),
                                    data=x)
            rs=cor(predict(rpi,newdata=x),x[,i],use="pairwise.complete.obs",method="spearman")^2
            if (rs < sum(rss)+rs.threshold) {
                break
            }
            A[i,vimps[vi]]=rs-sum(rss)
            rss=c(rss,rs-sum(rss))
        }
    }
    RS=A
    A[A>0]=1
    #A[C<0 & A > 1]=-1
    A=A+t(A)
    A[A>1]=1
    A[A==1 & C < 0]=-1
    if (keep=="mcgraph") {
        G=mcg.new(A,type=type)
        attr(G,"mode")="undirected"
        attr(G,"type")=type
        attr(G,"r.squared")=RS
        return(G)
    } else if (keep == "matrix") {
        return(A)
    } else if (keep == "list") {
        l=list()
        for (i in 1:ncol(A)) {
            l[[colnames(A)[i]]]=colnames(A)[which(A[i,] != 0)]
        }
        return(l)
    }
}

# https://stats.stackexchange.com/questions/14853/variable-importance-from-glmnet/211396#211396
# Very) Long story short, I advise to use the Agresti method:

# if X is the input matrix of the glmnet function,
# and cv.result is your glmnet object:
# sds <- apply(X, 2, sd)
# cs <- as.matrix(coef(cv.result, s = "lambda.min"))
# std_coefs <- coefs[-1, 1] * sds

# https://stats.stackexchange.com/questions/266592/how-to-calculate-r2-for-lasso-glmnet

#' @title Using Ridge, Elastic Net or Lasso regression for graph construction
#' @description The function `mcg.glmnet` - creates a new mcgraph object from
#'        given data using Ridge, Elastic Net or Lasso regressions. Nodes will be connected to other nodes if their coefficient is not zero or their R-square value
#'        exceeds a certain threshold.
#' @param x data frame or matrix with numerical daattr(Gswiss,"r.squared"),3ta.
#' @param type custom type name for the graph, Default: 'lasso'
#' @param rs threshold for the R-square value. Default: 0.04
#' @param alpha elasticnet mixing paramater with 0 we have Ridge regression, with 1 with have Lasso, in between we have elasticnet, default: 1
#' @return an mcgraph object for the data, or a adjacency matrix or list depending on the keep argument,
#'      in the adjacency matrix a value of minus indicate negative associations
#'      between the nodes, a value of +1 positive associations, additional attributes such as 'r.squared' values and 'std.coef'
#'      standardized coefficients are added as well
#' @examples
#' data(swiss)
#' Gswiss=mcg.glmnet(swiss,rs=0.1)
#' plot(Gswiss,main="Predicted Graph",layout='star')
# @importFrom glmnet glmnet cv.glmnet
#' @rdname mcg.glmnet
#' @author Detlef Groth, University of Potsdam
#' @export

mcg.glmnet <- function (x,type="lasso",rs=0.04,alpha=1) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package \"glment\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if (is.data.frame(x)) {
        x=as.matrix(x)
    }
    A=matrix(0,nrow=ncol(x),ncol=ncol(x))
    rownames(A)=colnames(A)=colnames(x)
    Coef=A
    rsqs=c()
    for (i in 1:ncol(x)) {
        xi=x[,-i]
        yi=x[,i]
        CV=glmnet::cv.glmnet(as.matrix(xi),yi,alpha=alpha,family="gaussian")
        lasso.model <- glmnet::glmnet(as.matrix(xi),yi, family = "gaussian",
                                      alpha=alpha,  lambda = CV$lambda.1se )
        rsqs=c(rsqs,lasso.model$dev.ratio)
        # normalize coefficients Agresti method(?)
        sds=apply(as.matrix(xi),2,sd)
        cs = as.matrix(coef(lasso.model, s = "lambda.1se"))
        std.coef = cs[-1, 1] * sds
        Coef[i,names(std.coef)]=std.coef
        c=coef(CV)[,1]
        c=c[2:length(c)]
        nm=names(which(c!=0))
        A[i,nm]=c[nm]
    }
    names(rsqs)=colnames(x)
    rsq=(abs(Coef)/apply(abs(Coef),1,sum))*rsqs
    C=A
    A[A>0]=1
    A[A<0]=-1
    A[rsq<rs]=0
    if (alpha == 1) {
        type="lasso"
    } else if (alpha == 0) {
        type="ridge"
    } else {
        type="elnet"
    }
    G=mcg.new(A,type=type)
    attr(G,"mode")="undirected"
    attr(G,"type")=type
    attr(G,"coef")=C
    attr(G,"std.coef")=Coef
    attr(G,"r.squared")=rsq
    return(G)
}

#' @title Determining the prediction quality by comparing predicted and true graph.
#'
#' @description  The function `mcg.accuracy` measures the prediction
#'  quality of a predicted graph in comparison to a known true one.
#'  The graph comparisons are done on the basis of undirected graphs only.
#'  Directed graphs are converted to undirected internally.
#'
#' @param g.true mcg.graph object or adjacency matrix of a the true graph.
#' @param g.pred mcg.graph object or adjacency matrix of a predicted graph.
#' @return returns list object with various accuracy measures, such as Sens(itivity), Spec(ificity), Acc(uracy), balanced classification rate (BCR), F1 measure and Mathhews correlaition coefficient (MCC).
#' @author: Detlef Groth <dgroth@uni-potsdam.de>
#' @examples
#' ang=mcg.angie(nodes=12,edges=18)
#' data=mcg.graph2data(ang)
#' pred=mcg.ct(t(data),rs=0.1)
#' round(unlist(mcg.accuracy(ang,pred)),2)
#' pred=mcg.lvs(t(data))
#' round(unlist(mcg.accuracy(ang,pred)),2)
#' pred=mcg.rpart(t(data))
#' round(unlist(mcg.accuracy(ang,pred)),2)
#' @export

mcg.accuracy <- function (g.true,g.pred) {
    # Convert to undirected
    g.true <- mcg.d2u(g.true)
    g.pred <- mcg.d2u(g.pred)
    g.true <- g.true[upper.tri(g.true)]
    g.pred <- g.pred[upper.tri(g.pred)]
    # Convert explictly to 'double' to avoid integer overflow
    TP=as.double(length(which(g.true != 0 & g.pred != 0)))
    FP=as.double(length(which(g.true == 0 & g.pred != 0)))
    TN=as.double(length(which(g.true == 0 & g.pred == 0)))
    FN=as.double(length(which(g.true != 0 & g.pred == 0)))
    # Inlcuded NA check
    Acc=if ( (TP+TN+FP+FN) != 0 ) { TP+TN/(TP+TN+FP+FN) } else { NaN }
    Sens=if ( (TP+FN) != 0 ) { TP/(TP+FN) } else { NaN }
    Spec=if ( (TN+FP) != 0 ) { TN/(TN+FP) } else { NaN }
    BCR=(Sens+Spec)/2
    F1=if ( (2*TP+FP+FN) != 0 ) { 2*TP/(2*TP+FP+FN) } else { NaN }
    MCC=if ( sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) != 0) { (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) } else { NaN }
    norm_MCC=(MCC+1)/2
    return(list(TP=TP,FP=FP,TN=TN,FN=FN,Sens=Sens,Spec=Spec, BCR=BCR,F1=F1,MCC=MCC, norm_MCC=norm_MCC))
}

#' @title Visualizing graph objects
#'
#' @description Plot function to visualize \code{mcgraph} objects.
#'
#' @details This is the basic plot functionality of \code{mcgraph} for plotting
#'          undirected and directed graphs. It is possible to give a layout
#'          for the graph by specifying the \code{layout} argument by a string
#'          or a two column matrix with the \code{x} and \code{y} coordinates
#'          for each node. The string must be one of the possible layout
#'          types given in \code{mcg.layout}, i.e. \code{mds} (multi-dimensional scaling), \code{sam}, \code{frg}
#'          (Fruchterman-Reingold force-directed layout), \code{circle} and \code{grid} (see \code{?mcg.layout} for more details).
#'          By default, the \code{frg} layout is used.
#'
#'          If no \code{dimnames} are given, the nodes will be named by \code{mcg.autonames}.
#'          It is also possible to omit node labels completely by setting the \code{show.labels}
#'          to \code{FALSE} or to give costum \code{labels} in form of a string vector.
#'
#'          Next to specifying basic properties like the color and the size for nodes, edges and the labels,
#'          there is also an \code{interactive} mode for modifying the layout by hand and saving the layout.
#'          Label properties can be specified via corresponding arguments used in the \code{text} function of the traditional \code{graphics}.
#'
#'          The function offers basic plotting functionalities. For more sophisticated use-cases, the \code{igraph}
#'          package by G. Csardi might be more appropriate.
#'
#' @param x Object of \code{mcgraph} or matrix class
#  @param \dots Parameters which are forwarded to the plot.igraph if igraph was loaded (TODO)
#' @param layout Either a character string defining how to plot the graph with
#'               possible values are \code{"mds"}, \code{"sam"}, \code{"frg"}, \code{"circle"}
#'               and \code{"grid"} or a matrix with x and y coordinates for each node. Default: \code{"frg"}
#' @param noise Should be noise added to the layout? Sometimes useful if nodes are too close. Default: \code{FALSE}
#' @param vertex.size Size of nodes. If \code{vertex.symbol="rectangle"}, width of rectangles.  Default: \code{1}
#' @param vertex.length Height of nodes, if \code{vertex.symbol="rectangle"}. Default: \code{0.5}
#' @param vertex.color Color of nodes. Default: \code{"grey80"}
#' @param vertex.in.color Color of in-nodes in directed graphs. Default: \code{"light blue"}
#' @param vertex.out.color Color of out-nodes in directed graphs. Default: \code{"salmon"}
#' @param vertex.border.color Color of the border of nodes. Default: \code{"black"}
#' @param vertex.symbol Shape of node symbol. Possible values: \code{"circle"}, \code{"rectangle"}. Default: \code{"circle"}
#' @param edge.width Thickness of edges. Default: \code{1.5}
#' @param edge.color Color of positive edges. Default: "\code{grey50}"
#' @param edge.neg.color Color of negative edges. Default: \code{"red"}
#' @param edge.lty Linetype of edges. Default: \code{"solid"}
#' @param arr.length Length of arrows. Default: \code{0.08}
#' @param arr.angle Angle of arrows. Default: \code{12}
#' @param arr.mode Specifies where the arrow heads are plotted. Either \code{"midpoint"} or \code{"aligned"}. Default: \code{"aligned"}
#' @param show.labels Should node labels be printed? Default: \code{TRUE}
#' @param label.size Size of the node text labels. Default: \code{1}
#' @param label.color Color of label. Default: \code{NULL}
#' @param label.font Font type to be used. \code{1} = normal, \code{2} = bold, \code{3} = italic, \code{4} = bold-italic. Default: \code{NULL}
#' @param label.adj Justification of labels. \code{O} for left/bottom, \code{1} for rigth/top, \code{0.5} for centered. Default: \code{NULL}
#' @param label.pos Position specifier. \code{1} = below, \code{2} = left, \code{3} = above, \code{1} = right. Default: \code{NULL}
#' @param label.offset Distance of label from specified coordinate if \code{label.pos} is specified. Default: \code{0.5}
#' @param label.vfont Null for default font family or character vector of length 2 for Hershey vector fonts. Default: \code{NULL}
#' @param label.text Use costum labels instead of automatic ones. Default: \code{NULL}
#' @param bw Boolean Should only black and white color schemes used? Default: \code{FALSE}
#' @param interactive Switch into interactive mode where you can click in the graph
#'                    and move nodes with two clicks, first selecting the node, second
#'                    click gives thehe new coordinates for the node. Default: \code{FALSE}
#' @param ... Other graphical parameters.
#' @return Returns the layout of the plotted network (invisible).
#' @author: Detlef Groth <email: dgroth@uni-potsdam.de>, Masiar Novine <email: masiar.novine@gmail.com>
#' @import graphics
#' @examples
#' band=mcg.band(16)
#' mcb=mcg.graph(band,n=200,input='A1')
#' plot(mcb,layout='circle')
#'
#' # Larger network
#' mca <- mcg.angie(200, 320)
#' plot(mcb,layout='circle')
#' plot(mcb,layout='circle')
#' @export

# not in use: @param layout either a layout function or a matrix with x and y coordinates for each node

plot.mcgraph <- function(x, layout=NULL, noise=FALSE,
                        vertex.size=1, vertex.length=0.5, vertex.color="orange", vertex.in.color="light blue",
                        vertex.out.color="salmon", vertex.border.color="black", vertex.symbol="circle",
                        edge.width=1.5, edge.color="black", edge.neg.color="red", edge.lty="solid",
                        arr.length=0.08, arr.angle=12, arr.mode="aligned", show.labels=TRUE, label.size=1, label.font=NULL,
                        label.color=NULL, label.adj=NULL, label.pos=NULL, label.offset=0.5, label.text=NULL, label.vfont=NULL,
                        bw=FALSE, interactive=FALSE, ...) {
    theta <- x
    # Check for dimnames, mandatory, otherwise nothing works
    if (any(is.null(colnames(theta)), is.null(rownames(theta)))) {
        rownames(theta) <- mcg.autonames(LETTERS, nrow(theta))
    }
    # Calculate layout, if no layout is given
    if (is.null(layout) & is.null(attr(x, 'layout'))) {
        layout <- 'frg'
    # Is it a object with an attribute 'layout' to use, e.g. a 'mcgraph' object?
    } else if (is.null(layout) & !is.null(attr(x, 'layout'))) {
        layout <- attr(x, 'layout')
        colnames(layout) <- c("x", "y")
        rownames(layout) <- colnames(theta)
    }
    # Modified to also cut beginning of line
    marrow <- function(x0, y0, x1, y1, cut.front, cut.back, col='black', ...){
        x0.new <- (1 - cut.front) * x0 + cut.front * x1
        y0.new <- (1 - cut.front) * y0 + cut.front * y1
        x1.new <- (1 - cut.back) * x1 + cut.back * x0
        y1.new <- (1 - cut.back) * y1 + cut.back * y0
        arrows(x0.new, y0.new, x1.new, y1.new, col=col, ...)
    }
    # 'asp=1' is mandatory to get optimal aligned arrows
    doPlot <- function (theta, xy, node=NULL, asp=1, ...) {
        if (!all(t(theta) == theta)) {
            graph_mode <- "directed"
        } else {
            graph_mode <- "undirected"
        }
        if (bw) {
            vertex.color <- "grey90"
            if (graph_mode == "directed") {
                vertex.out.color <- "grey45"
                vertex.in.color <- "grey65"
            }
        }
        xrange <- range(xy[, 1])
        yrange <- range(xy[, 2])
        xlim <- c(xrange[1] - 1/10 * diff(xrange), xrange[2] + 1/10 * diff(xrange))
        ylim <- c(yrange[1] - 1/10 * diff(yrange), yrange[2] + 1/10 * diff(yrange))
        # Empty plot
        plot(xy, type="n", axes=FALSE, xlab="", ylab="", xlim=xlim, ylim=ylim, asp=asp, ...)
        # Vectorized part
        if (graph_mode == "directed") {
            idx <- which(theta != 0)
        } else {
            # Get only the upper.tri indices, this avoids plotting edges twice
            idx <- which(upper.tri(theta) & theta != 0)
        }
        trg <- ceiling(idx / nrow(theta))
        src <- idx - nrow(theta) * (trg - 1)
        # Init vector for edge colors & adapt colors for negative edges
        edge.color <- rep(edge.color, length(idx))
        idx_neg <- which(theta < 0)
        edge.color[idx %in% idx_neg] = edge.neg.color
        x0 <- xy[rownames(theta)[src], 1]
        x1 <- xy[rownames(theta)[trg], 1]
        y0 <- xy[rownames(theta)[src], 2]
        y1 <- xy[rownames(theta)[trg], 2]
        # Solved: Alignment of arrow heads
        if (graph_mode == "undirected") {
            arrows(x0=x0, y0=y0, x1=x1, y1=y1, length=0, lwd=edge.width, col=edge.color, lty=edge.lty)
        }
        if (vertex.symbol == "circle") {
            # Multiply vertex size by factor to set default 'vertex.size' to 1
            vertex.size <- vertex.size * 1/20
            vertex.color <- rep(vertex.color, nrow(x))
            ins <- which(degree(x, mode="in") != 0 & degree(x, mode="out") == 0)
            vertex.color[ins] <- vertex.in.color
            outs <- which(degree(x, mode="out") != 0 & degree(x, mode="in") == 0)
            vertex.color[outs] <- vertex.out.color
            if (graph_mode == "directed") {
                r <- vertex.size
                x.user <- c(x0, x1)
                y.user <- c(y0, y1)
                # Euclidean distance
                xy.euc <- sqrt(diff(x.user, lag=length(x0))^2 + diff(y.user, lag=length(x0))^2)
                # Optimal would be to plot arrow heads independent of edges, but this seems to be too much effort, alternatively include 'shapes' library
                if (arr.mode == "midpoint") {
                    marrow(x0, y0, x1, y1, cut.front=0.45, cut.back=r / xy.euc, code=2, length=arr.length, angle=arr.angle, lwd=edge.width, col=edge.color, lty=edge.lty)
                    marrow(x0, y0, x1, y1, cut.front=r / xy.euc, cut.back=0.55, length=0, angle=arr.angle, lwd=edge.width, col=edge.color, lty=edge.lty)
                } else if (arr.mode == "aligned") {
                    marrow(x0, y0, x1, y1, code=2, length=arr.length, angle=arr.angle, lwd=edge.width, cut.front=r / xy.euc, cut.back=r / xy.euc, col=edge.color, lty=edge.lty)
                }
            }
            # Symbols with 'inches=FALSE' automatically adapts vertex size when the device window is resized
            symbols(x=xy[,1], y=xy[,2], bg=vertex.color, fg=vertex.border.color, circles=rep(vertex.size, nrow(xy)), add=TRUE, inches=FALSE)
        # MN Solved: arrow heads for rectangles
        } else if (vertex.symbol == "rectangle") {
            # This is just considered for 'rectangle'
            vertex.size <- vertex.size * 1/4
            vertex.length <- vertex.length * 1/4
            rect_dim <- matrix(c(rep(vertex.size, nrow(xy)), rep(vertex.length, nrow(xy))), nrow=nrow(xy), ncol=2)
            if (graph_mode == "directed") {
                delta_x = x1 - x0
                delta_y = y1 - y0
                phi = atan2(delta_y, delta_x) * 180 / pi
                x0.new = rep(0, times=length(phi))
                x1.new = rep(0, times=length(phi))
                y0.new = rep(0, times=length(phi))
                y1.new = rep(0, times=length(phi))
                for (i in seq_along(phi)) {
                    if (phi[i] < 45 & phi[i] > -45) {
                        x0.new[i] = x0[i] + (vertex.size / 2)
                        x1.new[i] = x1[i] - (vertex.size / 2)
                        y0.new[i]  = y0[i]
                        y1.new[i]  = y1[i]
                    } else if (phi[i] > 135 | phi[i] < -135) {
                        x0.new[i] = x0[i] - (vertex.size / 2)
                        x1.new[i] = x1[i] + (vertex.size / 2)
                        y0.new[i]  = y0[i]
                        y1.new[i]  = y1[i]
                    } else if (phi[i] >= 45 & phi[i] <= 135) {
                        x0.new[i] = x0[i]
                        x1.new[i] = x1[i]
                        y0.new[i]  = y0[i] + (vertex.length / 2)
                        y1.new[i]  = y1[i] - (vertex.length / 2)
                    } else if (phi[i] >= -135 & phi[i] <= -45) {
                        x0.new[i] = x0[i]
                        x1.new[i] = x1[i]
                        y0.new[i]  = y0[i] - (vertex.length / 2)
                        y1.new[i]  = y1[i] + (vertex.length / 2)
                    }
                }
                if (arr.mode == "midpoint") {
                    marrow(x0.new, y0.new, x1.new, y1.new, cut.front=0, cut.back=0.45, code=2, length=arr.length, angle=arr.angle, lwd=edge.width, col=edge.color, lty=edge.lty)
                    marrow(x0.new, y0.new, x1.new, y1.new, cut.front=0.55, cut.back=0, length=0, angle=arr.angle, lwd=edge.width, col=edge.color, lty=edge.lty)
                } else if (arr.mode == "aligned") {
                    arrows(x0=x0.new, y0=y0.new, x1=x1.new, y1=y1.new, length=arr.length, lwd=edge.width, code=2, angle=arr.angle)
                }
            }
            symbols(x=xy[,1], y=xy[,2], bg=vertex.color, fg=vertex.border.color, rectangles=rect_dim, add=TRUE, inches=FALSE)
        }
        # MN Costumized labels
        if (show.labels) {
            if (!is.null(label.text)) {
                text(x=xy, labels=label.text, adj=label.adj, pos=label.pos, offset=label.offset, vfont=label.vfont, cex=label.size, col=label.color, font=label.font, ...)
            }
            else {
                text(x=xy, labels=rownames(theta), adj=label.adj, pos=label.pos, offset=label.offset, vfont=label.vfont, cex=label.size, col=label.color, font=label.font, ...)
            }
        }
        invisible(xy)
    }
    if (any(class(layout) %in% "character")) {
       xy <- mcg.layout(x, mode=layout, noise=noise)
    } else {
        if (!any(class(layout) %in% "matrix")) {
            stop("given layout must be a character string such as mds or circle or a numeric matrix with two columns for x and y coordinates for each node")
        }
        if (ncol(layout) != 2) {
            stop("given layout matrix must have two columns")
        }
        if (nrow(layout) != nrow(theta)) {
            stop("given layout matrix must have same number of rows as the graph has nodes")
        }
        xy <- layout
    }
    # MN: Layout coordinates should be normalized in the range of -1 to 1, so that 'vertex.size' behaves predictable
    xy <- apply(xy, 2, function(x) -1 + ((x - min(x)) * (1 + 1)) / (max(x) - min(x)) )
    doPlot(theta, xy, ...)
    lay <- xy
    if (interactive) {
        print("Click two times: First on the point to move, second where to move. End with two right clicks!")
        while (TRUE) {
            loc <- locator(2)
            if (class(loc) == "NULL" | class(loc$x[2]) == "NULL" | class(loc$x[1]) == "NULL") {
                break
            }
            dlay <- rbind(lay, c(loc$x[1], loc$y[1]))
            d <- as.matrix(dist(dlay))[nrow(dlay), 1:(nrow(dlay) - 1)]
            nm <- names(which(d == min(d)))[1]
            lay[nm, 1] <- loc$x[2]
            lay[nm, 2] <- loc$y[2]
            lay <- doPlot(theta, lay, ...)
        }
    }
    invisible(lay)
}

#' @title summary.mcgraph gives back summary of mcgraph properties
#'
#' @description `summary.mcgraph` gives back summary of mcgraph properties.
#'
#' @details  This function gives back the graph properties such us type of graph, number of edges and nodes.
#' @param object object of class mcgraph.
#' @param \dots parameters which are forwarded to the generic summary functions. Currently not used.
#' @author: Detlef Groth <dgroth@uni-potsdam.de>
#' @examples
#' band=mcg.band(16)
#' summary(band)
#' @export

summary.mcgraph = function (object, ...) {
    x=object
    if (attr(x,"mode")=="undirected") {
        ecount=sum(x)/2
        density=ecount/(((nrow(x)^2)-nrow(x))/2)
    } else {
        ecount=sum(x)
        density=ecount/((nrow(x)^2)-nrow(x))
    }
    return(list(class=class(x),
                type=attr(x,"type"),
                mode=attr(x,'mode'),
                vcount=nrow(x),
                ecount=ecount,
                density=density))
}
#' @title Graph density of `mcgraph` objects
#'
#' @description `density.mcgraph` gives back the graph density of the graph object.
#'
#' @details  The graph density of a graph is the ratio of the number of edges divided by the  number of possible edges.
#' @param x object of class mcgraph.
#' @param \dots parameters which are forwarded to the generic density function. Currently not used.
#' @author: Detlef Groth <dgroth@uni-potsdam.de>
#' @examples
#' bara=mcg.barabasi(nodes=16,m=2)
#' density(bara)
#' @export

density.mcgraph = function (x,...) {
    su=summary(x)
    return((su$ecount/((su$vcount^2 - su$vcount)/2)))
}

# install degree as S3 function
# without diasabling a possible pkg::degree function
if (try(class(degree),silent=TRUE)=="function") {
    if (!isS3stdGeneric("degree")) {
        degree.default=degree
    }
} else {
    degree.default=function(x,...) {
        pos=paste("package",class(x),sep=":")
        fun=try(get("degree",pos=pos),silent=TRUE)
        if (is.function(fun)) {
            fun(x,...)
        } else {
            stop(paste("Unkown class",class(x),"for function degree"))
        }
    }
}
#' @title Degree of `mcgraph` graphs
#'
#' @description Get the number of edges for vertices of a graph. Generic function.plots mcgraph objects
#'
#' @param x Object of class igraph, mcgraph or something else
#' @param \dots Arguments forwarded to the right namespace method.
#' @seealso \link[mcgraph:degree.mcgraph]{degree.mcgraph}
#' @export degree

degree = function (x,...) UseMethod("degree")

#' @title degree of the vertices
#' @description The degree of the vertex is the number of edges connecting the vertex to other vertices.
#' @param x The mcgraph object to analyze.
#' @param mode Character string. Ignored if the graph is undirected. For directed graphs "out" will give the outgoing edges, "in" will give the incoming edges and "all" will give the sum of in- and outgoing edges.
#' @param \dots Currently not used.
#' @return Named vector with the degree for each edge.
#' @details The degree of a vertex represents a fundamental property of a vertex. Vertices with degree 0 are unconnected to other vertices. The degree of a vertex is an important centrality measure.
#' @examples
#' ang=mcg.angie(nodes=12,edges=16)
#' anu=mcg.u2d(ang,input='A1')
#' degree(ang)
#' degree(anu)
#' degree(anu,mode="in")
#' @import MASS
#' @author Detlef Groth <email: dgroth@uni-potsdam.de>
#' @rdname degree.mcgraph
#' @export

degree.mcgraph <- function (x,mode='all',...) {
    A=abs(x)
    if (attr(A,"mode") == "undirected") {
        # Use colSums for transpose, way faster. M.N.
        #return(apply(A,1,sum))
        return(colSums(t(A)))
    } else if (mode == "all") {
        #return(apply(A,1,sum)+apply(A,2,sum))
        return(colSums(t(A)) + colSums(A))
    } else if (mode == "in") {
        #return(apply(A,2,sum))
        return(colSums(A))
    } else if (mode == "out") {
        #return(apply(A,1,sum))
        return(colSums(t(A)))
    }
}

#' @title Layout a graph for plotting
#' @description Layout a graph for plotting
#' @param A Adjacency matrix
#' @param mode Character string for the plotting algorithm. Available are 'mds', 'sam', 'frg', 'circle' and 'grid', Default: sam
#' @param noise Boolean. Should noise be added to the coordinates. Useful if vertices overlap.
#' @param ... Other parameters concerning individual layout functions. See 'mcg.layout.frg'.
#' @return matrix with two columns for x and y coordinates.
#' @examples
#' ang=mcg.angie(nodes=12,edges=16)
#' xy=mcg.layout(ang,mode='mds')
#' plot(ang,layout=xy)
#' @rdname mcg.layout
#' @author Detlef Groth <email: dgroth@uni-potsdam.de>
#' @export

mcg.layout <- function (A, mode='sam', noise=FALSE, ...) {
    if (mode %in% c('mds','sam')) {
        if(attr(A,"mode") == "directed") {
            U=A+t(A)
            U[U!=0]=1
            A= U
        }
        A=connectComponents(A=A)
        sp=mcg.shortest.paths(A)
        xy=cmdscale(sp)
        #sp[]=sp+rnorm(length(sp),mean=0,sd=0.3)
        if (mode=='mds') {
            dxy=base::as.matrix(dist(xy))
            diag(dxy)=1
            idx=which(dxy<0.05,arr.ind=TRUE)
            if (nrow(idx)>1) {
                for (i in 1:nrow(idx)) {
                    # add noise to nodes on the same point
                    #print(i)
                    n=idx[i,1]
                    xy[n,1]=xy[n,1]+rnorm(1,mean=0,sd=0.1)
                    xy[n,2]=xy[n,2]+rnorm(1,mean=0,sd=0.1)
                }
            }
            return(xy)
        } else {
            xy=xy+jitter(xy)
            #xy=mcg.sam(sp)
            xy=MASS::sammon(sp,y=xy,trace=FALSE)$points
            # add some jitter
            #sp[]=sp+rnorm(nrow(projection)*2,sd=sd(projection)*0.1)
            #projection[]=
            #xy = MASS::sammon(sp)$points
        }
    } else if (mode == 'circle') {
        x=0
        y=0
        a=0.5
        b=0.5
        rad2deg <- function(rad) {(rad * 180) / (pi)}
        deg2rad <- function(deg) {(deg * pi) / (180)}
        nodes=rownames(A)
        xy=matrix(0,ncol=2,nrow=length(nodes))
        rownames(xy)=nodes
        for (i in 1:length(nodes)) {
            t=deg2rad((360/length(nodes))*(i-1))
            xp = a*cos(t)*0.75 + x;
            yp = b*sin(t)*0.75 + y;
            xy[nodes[i],]=c(xp,yp)
        }
    } else if (mode == 'star') {
        x=0
        y=0
        a=0.5
        b=0.5
        rad2deg <- function(rad) {(rad * 180) / (pi)}
        deg2rad <- function(deg) {(deg * pi) / (180)}
        nodes=rownames(A)
        xy=matrix(0,ncol=2,nrow=length(nodes))
        rownames(xy)=nodes
        xy[1,]=c(0.0,0.0)
        for (i in 2:length(nodes)) {
            t=deg2rad((360/(length(nodes)-1))*(i-2))
            xp = a*cos(t)*0.75 + x;
            yp = b*sin(t)*0.75 + y;
            xy[nodes[i],]=c(xp,yp)
        }
    } else if (mode == 'grid') {
        xy <- mcg.layout.grid(A, ...)
    } else if (mode == 'frg') {
        # xy = mcg.layout.frg(A)
        # MN Included C++ implementation
        xy = mcg.layout.frg(A, ...)
    } else {
        stop("unknown layout. Use mds, circle, grid or frg as layout")
    }
    xy=scale(xy)
    if (noise) {
        xy=xy+rnorm(length(xy),mean=0,sd=0.1)
    }
    # Naming dimnames
    colnames(xy) <- c("x", "y")
    if (any(is.null(colnames(A)), is.null(rownames(A)))) {
        rownames(xy) <- mcg.autonames(LETTERS, nrow(A))
    } else {
        rownames(xy) <- rownames(A)
    }
    return(xy)
}

#' @title Grid layout
#' @description Layout nodes on a grid
#' @param A Adjacency matrix
#' @param noise Boolean. Should noise be added to the coordinates? Default: FALSE
#' @param layout_matrix Optional matrix indicating positions of variables set as strings. Default: NULL
#' @details By default, nodes are automatically positioned on the grid. For defining node positions manually, it is possible to provide
#'          a layout matrix with positions of the variables. The variable names in the fields of the matrix must correspond to the
#'          \code{dimnames} of the adjacency matrix. Empty space between the variables can be introduced by setting the according field to \code{NA}.
#' @return Two column matrix for x and y coordinates.
#' @examples
#' # Network estimation
#' pred <- mcg.lvs(swiss)
#' # Automatically plotting nodes on a grid
#' plot(pred, vertex.symbol="rectangle", layout="grid", vertex.size=2.5, vertex.length=1.25)
#'
#' # Define costume positional matrix
#' lmat <- rbind(c(NA, "Catholic", NA),
#'               c("Fertility", NA, "Education"),
#'               c("Infant.Mortality", NA, "Agriculture"),
#'               c(NA, "Examination", NA))
#' # Plotting based on costume grid layout
#' lay <- mcg.layout.grid(pred, layout_matrix=lmat)
#' plot(pred, vertex.symbol="rectangle", vertex.size=2.5, vertex.length=1.25, layout=lay)
#' @rdname mcg.layout.grid
#' @author Detlef Groth <email: dgroth@uni-potsdam.de>, Masiar Novine <email: masiar.novine@gmail.com>
#' @export

mcg.layout.grid <- function(A, noise=FALSE, layout_matrix=NULL) {
    if (is.null(layout_matrix)) {
        n=nrow(A)
        xy=matrix(0,ncol=2,nrow=nrow(A))
        rownames(xy)=rownames(A)
        mody=ceiling(sqrt(n))
        x=0
        y=0
        for (r in rownames(A)) {
            if (x %% mody == 0) {
                y=y+1
                x=0
            }
            x=x+1
            xy[r,]=c(x,y)
        }
    # Order layout manually by giving relative positional matrix
    } else {
        var_names <- colnames(A)
        xy <- matrix(0, nrow=length(var_names), ncol=2)
        colnames(xy) <- c("x", "y")
        rownames(xy) <- var_names
        for (i in seq_along(var_names)) {
            # To get the right order for x, y and variables
            xy[i, 2:1] <- c(nrow(layout_matrix), ncol(layout_matrix)) - which(layout_matrix == var_names[i], arr.ind=TRUE, useNames=FALSE)
        }
    }
    return(xy)
}

#' @title Force-directed layouts for graphs
#' @description Generating graph layouts based on the Fruchterman-Reingold algorithm.
#' @param A Adjacency matrix.
#' @param iter Number of refinement iterations. Default: 500
#' @param W Width of the frame. Default: 10
#' @param L Length of the frame. Default: 10
#' @param temp_prop Proportion of the width of the frame used as temperature. Default: 0.2
#' @param force_prop Factor influencing the strength of repulsion and attraction. Higher values tend to seperate nodes more. Default: 1.5
#' @param quench_prop Initial proportion used as the basis in the cooling function. See 'Details' for more information. Default: 0.9
#' @param simmering_prop Proportion of the inital temperature used for simmering. Default: 0.05
#' @param seed Initial seed, can also be set to \code{NULL}. Default: 1234
#' @details The area of the frame is given by \code{W} * \code{L}, while the temperature t is given by \code{temp_prop} * \code{W}.
#'          If there is too much space between components of a graph, e.g. for small cluster networks, it can be useful to decrease the area
#'          by setting lower values for \code{W} and {L}, e.g. \code{W = 4} and {L = 4}.
#'
#'          Cooling is based on quenching and simmering as described in the original publication (see 'References'),
#'          i.e. rapid cooling of the temperature and then staying at a low constant temperature level to allow for minor optimizations.
#'          The lower constant is given by the initial temperature mulitplied by \code{simmering_prop}.
#'          Cooling is based on the following function:
#'              t[i] = t[i-1] * \code{quench_prop}^i
#'          where t[i] and t[i-1] are the temperature values at iteration i and i-1, respectively.
#'          The default values aim to give good layouts for different kind of graphs.
#'
#'          The constant 'k' described in the original publication
#'          can be influenced by \code{force_prop}, explicitly:
#'              k = \code{force_prop} * sqrt(area / n)
#'          where n is the number of nodes of the graph.
#'          Hence, by setting higher values for \code{force_prop} the repulsion of the nodes will be increased. For very high values (>> 100) the
#'          layout will converge to a circle layout.
#'
#'          The default number of iterations for refinement given by \code{iter} should be sufficient for most graphs to reach a global optimum.
#'          Sometimes, changing the \code{seed} can improve the overall layout.
#'
#'          The FRG algorithm is especially useful for symmetric graphs, such as graphs representing geometric objects,
#'          e.g. cubes, prisms etc. (see 'Examples').
#' @return Position matrix with two columns for 'x' and 'y' coordinates.
#' @references Fruchterman, T.M.J and Reingold, E.M., (1991), Graph Drawing by Force-directed Placement, Software–Practice and Experience, pp.1129-1164.
#' @examples
#' ang <- mcg.angie(nodes=25, edges=60)
#'
#' # Save layout for reuse, precalculate layout
#' xy <- mcg.layout(ang, mode='frg')
#' plot(ang, layout=xy, label.size=.7, vertex.size=.5)
#'
#' # Cluster example, layout created by plot function
#' C <- mcg.cluster(nodes=18, edges=24)
#' plot(C, layout="frg")
#'
#' # Cube example
#' cube <- matrix(c(0,1,0,1,0,0,1,0,
#'                  1,0,1,0,0,0,0,1,
#'                  0,1,0,1,0,1,0,0,
#'                  1,0,1,0,1,0,0,0,
#'                  0,0,0,1,0,1,1,0,
#'                  0,0,1,0,1,0,0,1,
#'                  1,0,0,0,1,0,0,1,
#'                  0,1,0,0,0,1,1,0), nrow=8, ncol=8)
#' cube <- mcg.new(cube)
#' rownames(cube) <- colnames(cube) <- mcg.autonames(LETTERS, nrow(cube))
#' plot(cube, layout="frg")
#'
#' # Pentagonal prism
#' penta <- matrix(c(0,1,1,0,0,0,0,0,1,0,
#'                   1,0,0,1,0,0,0,0,0,1,
#'                   1,0,0,1,1,0,0,0,0,0,
#'                   0,1,1,0,0,1,0,0,0,0,
#'                   0,0,1,0,0,1,1,0,0,0,
#'                   0,0,0,1,1,0,0,1,0,0,
#'                   0,0,0,0,1,0,0,1,1,0,
#'                   0,0,0,0,0,1,1,0,0,1,
#'                   1,0,0,0,0,0,1,0,0,1,
#'                   0,1,0,0,0,0,0,1,1,0), ncol=10, nrow=10)
#' penta <- mcg.new(penta)
#' rownames(penta) <- colnames(penta) <- mcg.autonames(LETTERS, nrow(penta))
#' plot(penta, layout="frg")
#' @rdname mcg.layout.frg
#' @author Masiar Novine <email: masiar.novine@gmail.com>, Detlef Groth <email: dgroth@uni-potsdam.de>
#' @export

# MN Included C++ implementation
# Wrapper for Fruchterman & Reingold
mcg.layout.frg <- function(A, iter=500, W=10, L=10, temp_prop=0.2, force_prop=1.5, quench_prop=0.9, simmering_prop=0.05, seed=1234) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    # MN: Call by reference, init matrix beforehand to avoid copies
    # Will be filled in C++ Armadillo
    pos <- matrix(runif(nrow(A) * 2, 0, 1), nrow=nrow(A), ncol=2)
    colnames(pos) <- c("x", "y")
    frgCpp(pos, A=A, iter=iter, W=W, L=L, temp_prop=temp_prop, force_prop=force_prop, quench_prop=quench_prop, simmering_prop=simmering_prop)
    # MN: Check, if at least one of the dimnames is empty, if not take existing dimname
    if (any(is.null(colnames(A)), is.null(rownames(A)))) {
        rownames(pos) <- mcg.autonames(LETTERS, nrow(A))
    } else {
        rownames(pos) <- rownames(A)
    }
    return(pos)
}

#' @title Type check for `mcgraph` objects
#' @description Checks, whether an object is of type `mcgraph`.
#' @param x Any R object.
#' @return A logical value. TRUE if argument is a mcgraph object.
#' @examples
#' band=mcg.band(nodes=12)
#' is.mcgraph(band)
#' is.mcgraph(1)
#' @rdname is.mcgraph
#' @author Detlef Groth <email: dgroth@uni-potsdam.de>
#' @export

is.mcgraph = function (x) {
    return(inherits(x, "mcgraph"))
}
#' @title Extract the adjacency matrix from a mcgraph graph object.
#' @description Extract the adjacency matrix from a mcgraph graph object and remove all other attributes.
#' @param x A mcgraph object.
#' @param \dots currently not used.
#' @return The adjacency matrix without other mcgraph attributes.
#' @examples
#' band=mcg.band(nodes=12)
#' class(band)
#' adj=as.matrix(band)
#' is.mcgraph(adj)
#' is.matrix(adj)
#' class(adj)
#' @rdname as.matrix.mcgraph
#' @author Detlef Groth <email: dgroth@uni-potsdam.de>
#' @export

as.matrix.mcgraph = function (x, ...) {
    M=x
    for (att in names(attributes(x))) {
        if (att %in% c("dim","dimnames")) {
            next
        }
        attributes(M)[[att]] = NULL
    }
    return(M)
}

connectComponents = function (A, use.null=FALSE, P=NULL) {
    A=as.matrix(A)
    A=A+t(A)
    A[A>0]=1
    # Speed up by avoiding double calculation, if already calculated by another function
    if (is.null(P)) { P <- mcg.shortest.paths(A) }
    if (!any(P==Inf)) {
        return(A)
    }
    comp=mcg.components(A=A, P=P)
    nodes=c()
    tab=table(comp)
    for (n in names(tab)) {
        c=names(which(comp==n))
        if (tab[[n]] > 2) {
            Am=A[c,c]
            # todo min
            deg=degree(mcg.new(Am))
            idx=which(deg>0)
            minval=min(deg[idx])
            idx=which(deg == minval)[1]
            node=c[idx]
        } else {
            node = c[1]
        }
        nodes=c(nodes,node)
    }
    if (use.null) {
        # not used yet
        A=cbind(A,NULL=rep(0,nrow(A)))
        A=rbind(A,NULL=rep(0,ncol(A)))
        A['NULL',nodes]=1
        A[nodes,'NULL']=1
    } else {
        A[nodes,nodes]=1
        diag(A)=0
    }
    return(A)

}

#' @title Get component ids from adjacency matrices or mcgraph objects
#' @description Checks the components of the given adjacency matrix or the given mcgraph objects and returns a unique numerical id for each components.
#' @param A An adjacency matrix or a mcgraph object.
#' @param P Matrix with the length of the shortest paths between each node.
#' @return vector with component ids and node names labeling the ids.
#' @examples
#' cls=mcg.cluster(nodes=12,edges=15,cluster=3)
#' class(cls)
#' plot(cls,layout='sam')
#' mcg.components(cls)
#' @rdname mcg.components
#' @author Detlef Groth <email: dgroth@uni-potsdam.de>
#' @export

mcg.components = function (A, P=NULL) {
    A=as.matrix(A)
    A=A+t(A)
    A[A>0]=1
    comp=c()
    # Speed up by avoiding double calculation, if already calculated by another function
    if (is.null(P)) { P <- mcg.shortest.paths(A) }
    nodes=rownames(A)
    x=1
    while (length(nodes) > 0) {
        n=nodes[1]
        idx=which(P[n,] < Inf)
        ncomp=rep(x,length(idx))
        names(ncomp)=rownames(P)[idx]
        comp=c(comp,ncomp)
        nodes=setdiff(nodes,rownames(P)[idx])
        x=x+1
    }
    return(comp[rownames(A)])
}

#' @title Get common metrics for a given network classifier, as well as the Receiver Operating Characteristic (ROC) or precision-recall (PR) curve.
#' @description \code{mcg.roc} can be used to evaluate the quality of network estimates.
#'              The function returns common evaluation metrics for given classification instances for a single network classifier and allows
#'              for plots of ROC and precision-recall curves.
#'
#' @details If the true network given by `G` is directed, it will be transformed to an undirected graph. Negative or weighted edges
#'          are not considered and are converted to 1-entries. The diagonal of the adjacency matrix is excluded. If the argument `plot='ROCC'` or
#'          `plot='PRC'` is set, the ROC curve or the precision-recall curve will be plotted, respectively.
#'          The argument 'fun' either takes a classifier function from the `mcgraph` package, e.g.:
#'          \describe{
#'          \item{\code{mcg.lvs}}{forward variable selection in linear models}
#'          \item{\code{mcg.rpart}}{variable selection using decision trees}
#'          \item{\code{mcg.ct}}{pruning based on hard correlation thresholding}
#'          }
#'          or a costum function either given by its function name or in the form of \code{function(x) \{ \}}. The given function
#'          will then be matched and called, after evaluating possible additional arguments given for \code{...}.
#'          There are three requirements for a given function
#'          \enumerate{
#'              \item{it must take a data object as a first argument with variables as columns and values as rows}
#'              \item{there must be an argument named `thr` for a classification threshold taking a single numeric value}
#'              \item{the given function must return an adjacency matrix}
#'          }
#'          Functions from other packages might be given in the form of wrappers to satisfy the above conditions.
#'
#'          \code{mcg.roc} returns a list containing the given classification thresholds \code{thr} and following metrics for each instance:
#'          The accuracy \code{acc}, sensitivity \code{sen}, specificity \code{spe}, precision \code{pre}, false-postive rate \code{fpr},
#'          Matthews Correlation Coefficient [MCC] \code{mcc}, the normed MCC \code{norm_mcc} calculated by \code{(MCC + 1) / 2} and the F1-score \code{f1}.
#'          Based on the trapezoid rule, the area under the ROC curve and the area under the precision-recall curve are returned.
#'          Instances, which cannot be calculated because they are mathematically not defined, are set to \code{NA} by convention.
#'
#'          Note, that the raw results are returned in order of the given classification instances. For internal calculation of the area under
#'          the curve and plotting, additional start- and endpoints are be included and the values are sorted in descending order
#'          after excluding \code{NA} entries. For convenient reuse for own plots, those values are also returned by the nested list
#'          entries \code{ROCC} and \code{PRC}, respectively.
#'
#' @param d Data frame with variables as columns and data values as rows.
#' @param G Adjacency matrix of the true network.
#' @param thresholds Vector of classification thresholds in increasing order.
#' @param plot Either 'ROCC' for the classical ROC curve (sensitivity vs. 1 - FPR) or 'PRC' for the precision-recall curve. Default: 'ROCC'
#' @param fun The name of the classifier function, either unquoted or with quotation marks. Default: mcg.lvs
#' @param ... additional arguments for 'fun'
#' @return a list containing values for several metrics (see Details).
#' @examples
#' # Build lion star sign network
#' L <- mcg.constellation(name="lion")
#'
#' # Create Monte-Carlo data for given network
#' set.seed(1234)
#' d <- as.data.frame(t(mcg.graph2data(L)))
#'
#' # Create & plot ROC curves using forward linear variable selection
#' mcg.roc(d=d, G=L, thresholds=seq(0, 1, length.out=11), fun=mcg.lvs)
#'
#' # Variable selection using decision tree
#' mcg.roc(d=d, G=L, thresholds=seq(0, 1, length.out=11), fun=mcg.rpart)
#'
#' # Hard correlation thresholding using R-square values
#' mcg.roc(d=d, G=L, thresholds=seq(0, 1, length.out=11), fun=mcg.ct)
#'
#' # Example for a costum function to perform hard correlation thresholding
#' hct <- function(d, thr=0.5) {
#'     # Init adjacency matrix
#'     A <- matrix(0, nrow=ncol(d), ncol=ncol(d))
#'     # Calculate correlation matrix
#'     cor_mt <- cor(d)
#'     # Exclude diagonal of correlation matrix
#'     diag(cor_mt) <- 0
#'     # Simple hard correlation thresholding
#'     A[which(cor_mt >= thr)] <- 1
#'     return(A)
#' }
#' # Get metrics and plot Precision-Recall curve
#' mcg.roc(d=d, G=L, thresholds=c(0.2, 0.4, 0.5), plot="PRC", fun=hct)
#'
#' @rdname mcg.roc
#' @author Masiar Novine <email: masiar.novine@gmail.com>
#' @export

mcg.roc <- function(d, G, thresholds, plot="ROCC", fun=mcg.lvs, ...) {
    fun <- match.fun(fun)
    thr <- "thr"
    if (!is.matrix(G)) {
        G <- as.matrix(G)
    }
    # Don't consider negative weight egdes
    G <- abs(G)
    # Convert directed to undirected graphs
    if (!all(t(G) == G)) {
        G <- G + t(G)
        G[which(G != 0)] <- 1
    }
    # Exclude diagonal from calculations
    G <- G[upper.tri(G)]
    # Make weighted graph binary
    G[G != 0] <- 1
    # Allocate list
    res <- list()
    nms.res <- c("thr", "tp", "tn", "fp", "fn", "acc", "sen", "spe", "pre", "fpr", "f1", "mcc", "mcc.norm", "bcr", "auc", "pr.auc")
    for (l in seq_along(nms.res)) {
        res[[l]] <- rep(0, length(thresholds))
    }
    names(res) <- nms.res
    # Note: Could cause trouble, if by chance some other argument is called 'rs'
    if (any(names(formals(fun)) == "rs")) { thr <- "rs" }
    for (i in seq_along(thresholds)) {
        formals(fun)[thr] <- thresholds[i]
        # Take only absolute values
        P <- abs(forceAndCall(1, fun, d, ...))
        # Just in case: Convert directed to undirected graphs
        if (!all(t(P) == P)) {
            P <- P + t(P)
            P[which(P != 0)] <- 1
        }
        # Make weighted graph binary
        P[P != 0] <- 1
        # Exclude diagonal of predictions
        P <- P[upper.tri(P)]
        # Explicitly convert to 'double' to avoid overflow regarding MCC
        tp <- as.double(length(intersect(which(G != 0), which(P != 0))))
        fp <- as.double(length(which(P != 0))) - tp
        tn <- as.double(length(intersect(which(G == 0), which(P == 0))))
        fn <- as.double(length(which(P == 0))) - tn
        res$thr[[i]] <- thresholds[i]
        # Check if denominator is not 0, else set to NaN to avoid Inf cases
        res$acc[[i]] <- if ( (tp + tn + fp + fn) != 0 ) { (tp + tn) / (tp + tn + fp + fn) } else { NaN }
        res$sen[[i]] <- if ( (tp + fn) != 0 ) { tp / (tp + fn) } else { NaN }
        res$spe[[i]] <- if ( (tn + fp) != 0 ) { tn / (tn + fp) } else { NaN }
        res$pre[[i]] <- if ( (tp + fp) != 0 ) { tp / (tp + fp) } else { NaN }
        res$f1[[i]] <- if ( (2 * tp + fp + fn) != 0 ) { (2 * tp) / (2 * tp + fp + fn) } else { NaN }
        res$mcc[[i]] <- if ( sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) != 0 ) { (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) } else { NaN }
        res$mcc.norm[[i]] <- (res$mcc[[i]] + 1) / 2
        res$bcr[[i]] <- 0.5 * res$sen[[i]] + 0.5 * res$spe[[i]]
        res$fpr[[i]] <- 1 - res$spe[[i]]
        res$tp[[i]] <- tp
        res$tn[[i]] <- tn
        res$fp[[i]] <- fp
        res$fn[[i]] <- fn
    }
    # Include additional first & last value for calculations & plotting
    # Exclude NAs
    sen_auc <- c(1, res$sen[which(!is.na(res$sen) & !is.na(res$fpr))], 0.0)
    fpr <- c(1, res$fpr[which(!is.na(res$sen) & !is.na(res$fpr))], 0.0)
    pre <- c(0.0, res$pre[which(!is.na(res$pre) & !is.na(res$sen))], 1)
    sen_pr_auc <- c(1, res$sen[which(!is.na(res$pre) & !is.na(res$sen))], 0.0)

    # AUC
    # Sort y values by x values to keep pair relationship
    names(fpr) <- seq_along(1:length(fpr))
    fpr <- sort(fpr, decreasing=TRUE)
    sen_auc <- sen_auc[strtoi(names(fpr))]

    # PR-AUC
    # Sort y values by x values to keep pair relationship
    names(sen_pr_auc) <- seq_along(1:length(sen_pr_auc))
    sen_pr_auc <- sort(sen_pr_auc, decreasing=TRUE)
    pre <- pre[strtoi(names(sen_pr_auc))]

    # Calculate AUC based on trapezoid rule
    delta_w1 <- sen_auc[1:(length(sen_auc) - 1)]
    delta_w2 <- sen_auc[2:(length(sen_auc))]
    h1 <- abs(diff(fpr))
    res$auc <- sum(((delta_w1 + delta_w2) / 2) * h1, na.rm=TRUE)

    # Calculate PR-AUC based on trapezoid rule
    delta_z1 <- pre[1:(length(pre) - 1)]
    delta_z2 <- pre[2:(length(pre))]
    h2 <-  abs(diff(sen_pr_auc))
    res$pr.auc <- sum(((delta_z1 + delta_z2) / 2) * h2, na.rm=TRUE)

    # Remove names
    names(fpr) <- names(sen_pr_auc) <- NULL

    # Include editted metrics for convenience
    res$ROCC$fpr <- fpr
    res$ROCC$sen <- sen_auc
    res$PRC$sen <- sen_pr_auc
    res$PRC$pre <- pre

    # For plotting
    if (!missing(plot)) {
        if (plot == "ROCC") {
            plot(x=fpr, y=sen_auc, type="l", lwd=1.5, pch=20, xlab="FPR", ylab="Sensitivity",
                panel.first=abline(h=seq(0, 1, length.out=11), v=seq(0, 1, length.out=11), col="lightgray", lty="dotted"))
            points(x=fpr[-c(1, length(fpr))], y=sen_auc[-c(1, length(sen_auc))], pch=21, bg="black")
            lines(x=seq(0, 1, length.out=5), y=seq(0,1,length.out=5), col="red", lty="dotted")
            text(0.5, 0.5, paste0("AUC = ", round(res$auc, 3)))
        }
        else if (plot == "PRC") {
            plot(x=sen_pr_auc, y=pre, type="l", lwd=1.5, pch=20, xlab="Sensitivity", ylab="Precision",
                panel.first=abline(h=seq(0, 1, length.out=11), v=seq(0, 1, length.out=11), col="lightgray", lty="dotted"))
            points(x=sen_pr_auc[-c(1, length(sen_pr_auc))], y=pre[-c(1, length(pre))], pch=21, bg="black")
            lines(x=seq(0, 1, length.out=5), y=seq(1, 0, length.out=5), col="red", lty="dotted")
            text(0.5, 0.5, paste0("PR-AUC = ", round(res$pr.auc, 3)))
        }
    }
    return(res)
}

#' @title Get the shortest paths in a graph
#' @description Returns the shortes paths between each pair of nodes or Inf if they belong
#'     to different graph components which are not connected. Attention this only works for undirected graphs and directed graphs are converted internally to undirected graphs.
#'     The method is mainly used for the layout mechanism.
#' @param A mcgraph object or adjacency matrix
#' @param mode should graph be taken as directed graphs, or undirected, if mode='undirected' is given, the graph is transformed into an undirected graph, if mode is 'directed' no graph transformation is done, default: 'directed'
#' @param code should the R or the C++ version be used, default: C++
#' @return matrix of same dimension as input matrix with vector with component ids and node names labeling the ids.
#' @examples
#' cls <- mcg.cluster(nodes=12, edges=14, cluster=3)
#' class(cls)
#' plot(cls, layout='sam')
#' mcg.shortest.paths(cls, code="R")
#' mcg.shortest.paths(cls, code="C++")
#' ang <- mcg.angie(node=36, edges=48)
#' t1 <- Sys.time()
#' mcg.shortest.paths(cls, code="R");
#' print(Sys.time() - t1)
#' t1 <- Sys.time()
#' mcg.shortest.paths(cls, code="C++")
#' print(Sys.time() - t1)
#' @rdname mcg.components
#' @author Detlef Groth <email: dgroth@uni-potsdam.de>, Masiar Novine <email: masiar.novine@gmail.com>
#' @export

mcg.shortest.paths = function (A,mode="directed", code="C++") {
    if (mode == "undirected") {
        A=A+t(A)
        A[A!=0]=1
    }
    if (code == "C++") {
        # Init output beforehand
        S <- matrix(Inf, nrow=nrow(A), ncol=ncol(A))
        diag(S) <- 0
        # MN call by reference
        bfsCpp(S, as.matrix(A))
        rownames(S)=colnames(S)=rownames(A)
        attr(S,"class")=attr(A,"class")
        attr(S,"mode")=attr(A,"mode")
        attr(S,"type")=attr(A,"type")
    } else if (code == "R") {
        S=A
        S[]=Inf
        diag(S)=0
        x=1
        S[A > 0 & A < Inf]=1
        while (TRUE) {
            flag = FALSE
            for (m in 1:nrow(S)) {
                ns=which(S[m,] == x)
                for (n in ns) {
                    for (o in which(A[n,]==1)) {
                        if (o != m) {
                            flag = TRUE
                            if (S[m,o] > x + 1) {
                                S[m,o]=x+1
                                if (mode == "undirected") {
                                    S[o,m]=x+1
                                }
                            }
                        }
                    }
                }
            }
            if (!flag) {
                break
            }
            x=x+1
        }
    } else {
        stop("Error: argument code must bei either 'C++' or 'R'!")
    }
    return(S)
}

#' @title Conversion of undirected to directed graphs
#'
#' @description `mcg.u2d` gets an undirected input graph and returns a directed graph based on one or more input nodes.
#'
#' @details This function creates a directed graph from an undirected one by the given input nodes. Input nodes can be chosen by names or a number for random selection of input nodes will be given.
#' Input nodes will have at least shortest path distance to other input nodes of pathlength two.
#' Selected input nodes will draw in each iteration outgoing edges to other nodes in the nth iteration neighborhood. The input
#' nodes will alternatively select the next edges on the path to not visited nodes. All edges will be only visited onces.
#'
#' @param A undirected input graph
#' @param input number of input nodes in the graph, if number of inout nodes is smaller than number of components, for each component one inout node is automatically created.
#' @param negative proportion of inhibitive associations in the network value between 0 and 1 are acceptable, Default 0.0
#' @return directed graph
#' @author: Detlef Groth <dgroth@uni-potsdam.de>
#' @examples
#' cross=mcg.cross(bands=5,length=5)
#' D=mcg.u2d(cross,input=3)
#' plot(D)
#' angie=mcg.angie(nodes=8,edges=12)
#' and=mcg.u2d(angie,input=2,negative=0.2)
#' plot(and,layout="sam")
#' hubs=mcg.hubs(nodes=18,hubs=3)
#' plot(mcg.u2d(hubs,input=4),layout="sam")
#' @export

mcg.u2d = function (A,input=2,negative=0.0) {
    # creates from undirected a directed one
    # A must be an symmetric adjacency matrix
    # MN The new check is 5x faster and takes less memory
    #if (!identical(A[lower.tri(A)], t(A)[lower.tri(A)])) {
    if (!all(A == t(A))) {
        stop("adjacency matrix must be symmetric")
    }
    if (negative < 0 || negative > 1) {
        stop("negative proportions must be within 0 and 1")
    }
    # undirected matrix
    U=A
    # future directed matrix
    D=U
    D[]=0
    neighbours=c()
    if (class(input)=='numeric') {
        nodes=c()

        comps=mcg.components(A)
        if (max(comps)>1) {
            for (i in 1:max(comps)) {
                if (length(which(comps==i)) > 1) {
                    node=sample(names(which(comps==i)),1)
                    neighbours=c(neighbours,rownames(U)[which(U[node,]==1)],node)
                    nodes=c(nodes,node)
                }
            }
        }
        n=input-length(nodes)
        while (n>0) {
            rnames=setdiff(rownames(U),c(neighbours))
            node=sample(rnames,1)
            n=n-1
            neighbours=c(neighbours,rownames(U)[which(U[node,]==1)],node)
            nodes=c(nodes,node)
        }
    } else {
        nodes=input
        n=0
    }
    visits=list()
    while (length(nodes)>0) {
        node=nodes[1]
        edges=which(U[node,]==1)
        newnodes=colnames(U)[edges]
        if (length(nodes)==1) {
            nodes=newnodes
        } else {
            nodes=c(nodes[2:length(nodes)],newnodes)
        }
        D[node,edges]=1
        U[node,]=0
        U[,node]=0
    }
    A[]=D
    if (negative>0) {
        idx=which(A==1)
        n=floor(length(idx)*negative)
        if(n>0) {
            min=sample(idx,n)
            A[min]=-1
        }
    }
    attr(A,"mode")="directed"
    return(A)
}

#' @title Conversion of directed to undirected graphs
#'
#' @description `mcg.d2u` gets a directed input graph and returns an undirected graph where every edges becomes undirected.
#'
#' @param A directed input graph
#' @return undirected mcgraph object
#' @author: Detlef Groth <dgroth@uni-potsdam.de>
#' @examples
#' cross=mcg.cross(bands=5,length=5)
#' D=mcg.u2d(cross,input=3)
#' U=mcg.d2u(D)
#' par(mfrow=c(1,2))
#' plot(D)
#' plot(U)
#' @export

mcg.d2u = function (A) {
    U=A+t(as.matrix(A))
    U[U>1]=1
    attr(U,'mode')="undirected"
    return(U)
}

#' @title Generation of Monte Carlo simulated data
#'
#' @description `mcg.graph2data` create simulated data from a directed graph
#'
#' @details This function takes a input either a adjacency matrix from a directed graph, or a directed mcgraph created with the mcg.u2d function.
#' The function creates simulated data for the given graph using a Monte Carlo simulation.
#' With the defaults the function will create data where the absolute correlations
#' are on average around 0.5, To increase or decrease those average correlations,
#' increase  or decrease the number of iterations using the 'iter' argument.
#'
#' @param A either a adjacency matrix from a directed graph or a mcgraph object created with `mcg.u2d`.
#' @param n number of samples for which data should be generated.
#' @param iter Number of iterations done for each node to generate the data. Default: 30
#' @param val initial value to be used to for nodes. Default: 100
#' @param sd standard deviation to be used at initialisation and in each iteration to create new values for each node. Default: 2
#' @param prop amount of influence of source node for input node per iteration, default: 0.05
#' @param noise amount of scatter added in each iteration to the computed value. It is the sd of a rnorm call with mean of zero. Default: 1
#' @param init initialization values, matrix should match the A adjacency matrix with rows and the requested number of samples as columns, if NULL start values are choosen randomly default: NULL
#' @param code should the R or the C++ version be used, default: C++
#' @return data matrix where rows are the node names and columns are the values for the simulated data. Data are not normalized to allow inspection of the data influence. Before further usage data should be scaled.
#'
#' @examples
#' anu <- mcg.angie(nodes=12, edges=18)
#' and <- mcg.u2d(anu, input=2, negative=0.2)
#' and.data <- mcg.graph2data(and, n=200)
#' t1 <- Sys.time()
#' and.data <- mcg.graph2data(and, n=200, code="R")
#' print(Sys.time() - t1)
#' t1 <- Sys.time()
#' and.data <- mcg.graph2data(and, n=200, code="C++")
#' print(Sys.time() - t1)
#' round(cor(t(and.data)), 2)
#' plot(and)
#' # directed graph example
#' G=matrix(0,nrow=7,ncol=7)
#' rownames(G)=colnames(G)=LETTERS[1:7]
#' G['A','C']=1
#' G['B','C']=1
#' G['C','D']=1
#' G['D','E']=1
#' G['D','F']=1
#' G['E','F']=1
#' G['G','A']=-1
#' G
#' g=mcg.new(G)
#' plot(g,layout='sam')
#' data=mcg.graph2data(g)
#' print(round(cor(data),1))
#' # weighted graph
#' G['A','C']=1.6
#' g=mcg.new(G)
#' plot(g,layout='sam')
#' data=mcg.graph2data(g)
#' print(round(cor(data),1))
#' @keywords simulation
#' @import stats
#' @author Detlef Groth <email: dgroth@uni-potsdam.de>, Masiar Novine <email: masiar.novine@gmail.com>
#' @export

mcg.graph2data <- function (A,n=200,iter=30,val=100,sd=2,prop=0.05,noise=1,init=NULL,code="C++") {
    if (code == "C++") {
        # DG: TODO (see vignette example) C++ does not work yet!!
        # MN Fixed issue of init: I included a boolean for the C++ version to differentiate between given init and generated matrix
        if (is.null(init)) {
            data <- matrix(0, nrow=ncol(A), ncol=n)
            bool_init <- FALSE
        } else {
            data <- init
            bool_init <- TRUE
        }
        # MN Avoid copy by using call by reference of the resulting matrix
        graph2dataCpp(D=data, A=A, n=n, iter=iter, val=val, sd=sd, prop=prop, noise=noise, init=bool_init)
        rownames(data) <- rownames(A)
        return(data)
    } else if (code != "R") {
        stop("Error: code must be either R or C++")
    }
    res=matrix(0,ncol=n,nrow=nrow(A))
    rownames(res)=rownames(A)
    #keep=list()
    for (i in 1:n) {
        if (class(init[1]) == "NULL") {
            units=rnorm(nrow(A),mean=val,sd=sd)
        } else {
            units=init[,i]
        }
        names(units)=rownames(A)
        for (j in 1:iter) {
            for (node in sample(rownames(A))) {
                targets=colnames(A)[which(A[node,]!=0)]
                for (target in sample(targets)) {
                    if (TRUE) {
                        P=abs(A[node,target])
                        nval=units[[node]]*(prop*P)
                        nval=nval+units[[target]]*(1-(prop*P))
                        if (A[node,target]<0) {
                            diff=nval-units[[target]]
                            nval=units[[target]]-diff
                        }
                        units[[target]]=nval;
                    } else {
                        units[[target]]=units[[target]]+(units[node]*prop*A[node,target])
                    }

                }
            }
            units=units+rnorm(length(units),sd=noise)
        }

        res[,i]=units
    }
    return(res)
}

#' @title Imputation of missing values using rpart, knn, mean or median
#' @description `mcg.impute` imputes missing data in data frames or matrices using either
#'   decision trees, k-nearest-neighbor approach with Eucildean distances or
#'   simple mean and median computations for the variables.
#' @section Details:
#'     Many mathematical methods did not allow missing data in their inputs.
#'     The missing values have to be guessed in this case.
#'     The more basic approaches are replacing the NA values with the median or mean for the variable.
#'     More advanced methods are using the mean only for a few samples which are very similar to the sample where the values is missing.
#'     This method is called knn-imputation and uses the k-nearest neighbors (default here 5) only for computing the replacement value.
#'     For data where the value order in the columns is of importance it is often desired to replace missing values with the mean of the two closest neighbors in the data frame or matrix, for instance if the data are ordered by time. Here the method timemean follows this approach.
#'     The method rpart uses decision trees to impute the values, it is currently the only method to impute as well factor variables.  The advantage of the rpart method is that it can impute not only numerical values but as well factor variables.
#'
#' @param data the data frame or matrix with missing values.
#' @param method which method to be used to impute the values, either 'mean', 'median', rpart, 'knn', default: 'rpart'
#' @param k how many nearest neighbors to be used if method is 'knn', default: 5
#' @return data with imputed values.
#' @examples
#' data(iris)
#' set.seed(123)
#' ir=as.matrix(iris[,1:4])
#' ir.mv=ir
#' # introduce 5 percent NA's
#' mv=sample(1:length(ir),as.integer(0.05*length(ir)))
#' ir.mv[mv]=NA
#' ir.imp.med=mcg.impute(ir.mv,method='median') # not good
#' ir.imp.rpart=mcg.impute(ir.mv) # method rpart (default)
#' ir.imp.knn=mcg.impute(ir.mv,method='knn')
#' rmse = function (x,y) { return(sqrt(sum((x-y)^2))) }
#' rmse(ir[mv],ir.imp.med[mv]) # should be high
#' rmse(ir[mv],ir.imp.rpart[mv]) # should be low!
#' rmse(ir[mv],ir.imp.knn[mv]) # should be low!
#' cor(ir[mv],ir.imp.med[mv])
#' cor(ir[mv],ir.imp.rpart[mv])
#' cor(ir[mv],ir.imp.knn[mv]) # should be high!
#' # factor variables
#' data(iris)
#' ciris=iris
#' idx=sample(1:nrow(ciris),15) # 10 percent NA's
#' ciris$Species[idx]=NA
#' summary(ciris)
#' ciris=mcg.impute(ciris,method="rpart")
#' table(ciris$Species[idx],iris$Species[idx])
#' @name mcg.impute
#' @import rpart
#' @export
#' @seealso \link{mcg.new}.
#' @author Detlef Groth <email: dgroth@uni-potsdam.de>

mcg.impute <- function (data,method='rpart',k=5) {
    if (method %in% c('mean','median')) {
        if (is.data.frame(data) | is.matrix(data)) {
            for (i in 1:ncol(data)) {
                data[,i]=mcg.impute(data[,i],method=method)
            }
        } else {
            if (method == "median") {
                mn=median(data,na.rm=TRUE)
            } else {
                mn=mean(data,na.rm=TRUE)
            }
            idx=which(is.na(data))
            data[idx]=mn
        }
        return(data)
    } else if (method == "rpart") {
        # TODO: refinement for many variables,
        # take only variables with highest absolute correlation
        # into account if more than 10 variables take top 10
        idata=data
        for (i in 1:ncol(data)) {
            idx = which(!is.na(data[,i]))
            if (length(idx) == nrow(data)) {
                next
            }
            if (is.factor(data[,i])) {
                model=rpart::rpart(formula(paste(colnames(data)[i],"~.")),
                            data=as.data.frame(data[idx,]),
                            method="class")
                x2 = predict(model,newdata=as.data.frame(data[-idx,]),
                             type="class")
            } else {
                model=rpart::rpart(formula(paste(colnames(data)[i],"~.")),
                            data=as.data.frame(data[idx,]))
                x2 = predict(model,newdata=as.data.frame(data[-idx,]))
            }

            idata[-idx,i]=x2
        }
        return(idata)
    } else if (method == "knn") {
        if (ncol(data) < 4) {
            stop("knn needs at least 4 variables / columns")
        }
        data.imp=data
        #D=as.matrix(1-((cor(t(data.imp),use="pairwise.complete.obs")+1)/2))
        D=as.matrix(dist(scale(data.imp)))
        for (i in 1:ncol(data)) {
            idx=which(is.na(data[,i]))
            idxd=which(!is.na(data[,i]))
            for (j in idx) {
                idxo=order(D[j,])
                idxo=intersect(idxo,idxd)
                mn=mean(data[idxo[1:k],i])
                data.imp[j,i]=mn
            }
        }
        return(data.imp)
    } else {
        stop("Unkown method, valid methods are rpart, knn, mean or median")
    }
}

#' @title Creation of band or chain graph
#'
#' @description `mcg.band` creates a band, alias chain graph with given number of nodes
#'
#' @details Create a simple band or chain graph with a given number of nodes.
#' @param nodes number of nodes to create.
#' @return graph of class mcgraph with adjacency matrix for an undirected graph
#' @author Detlef Groth <dgroth@uni-potsdam.de>
#' @examples
#' band=mcg.band(12)
#' plot(band,vertex.color="salmon")
#' @export

mcg.band = function (nodes=8) {
    mt=matrix(0,nrow=nodes,ncol=nodes)
    rownames(mt)=colnames(mt)=mcg.autonames(LETTERS,nodes)
    for (row in 1:(nrow(mt)-1)) {
        mt[row,row+1]=1
        mt[row+1,row]=1
    }
    class(mt)="mcgraph"
    attr(mt,"type")="band"
    attr(mt,"mode")="undirected"
    return(mt)
}

#' @title Example networks based on star sign constellations
#'
#' @description `mcg.constellation` creates graphs, mimicking selected start sign constellations.
#'
#' @details This graph generation function generates graphs which
#'          mimick constellations on the northern sky, e. g. 'Lion',
#'          'Dipper' and 'Virgin'. There will be a few default input nodes given.
#' @param name the name of the constellation, i. e. 'Lion', Dipper' or 'Virgin'. Default: 'Lion'
#' @return graph of class mcgraph with adjacency matrix for an undirected graph
#' @author Masiar Novine <email: masiar.novine@gmail.com>
#' @examples
#' lion <- mcg.constellation(name='lion')
#' # Plot with adjusted vertex, edge and label size for nicer look
#' plot(lion, main="Lion", label.size=0.7, vertex.size=0.8)
#' @export

mcg.constellation <- function (name="lion") {
    name=tolower(name)
    if (name == "lion") {
        nodes <- 13
        A <- matrix(0, nrow=nodes, ncol=nodes)
        rownames(A) <- colnames(A) <- mcg.autonames(LETTERS, nodes)
        st <- c('A1', 'H1', 'H1', 'E1', 'L1', 'K1', 'M1', 'M1', 'F1','G1', 'G1', "D1", "D1", "B1", "I1", "J1")
        ed <- c('H1', 'I1', 'E1', 'L1', 'K1', 'M1', 'E1', 'F1', 'G1','H1', 'D1', "I1", "B1", "I1", "J1", "C1")
        for (i in 1:length(st)) {
            A[st[i], ed[i]] <- 1
        }
        leo_lay <- matrix(c(
             0.48, 0.21,
            -0.65, 0.23,
            -0.13, 0.17,
            -0.16, 0.27,
             0.69, 0.29,
             0.37, 0.30,
             0.31, 0.27,
             0.49, 0.25,
            -0.16, 0.24,
            -0.22, 0.21,
             1.00, 0.32,
             0.95, 0.29,
             0.61, 0.31),
            nrow=13, byrow=TRUE)
        rownames(leo_lay) <- rownames(A)
        colnames(leo_lay) <- c("x", "y")
        class(A) <- "mcgraph"
        attr(A,"type") <- "constellation"
        attr(A,"mode") <- "directed"
        attr(A,"layout") <- leo_lay
    } else if (name == "dipper") {
        nodes <- 7
        A <- matrix(0, nrow=nodes, ncol=nodes)
        rownames(A) <- colnames(A) <- mcg.autonames(LETTERS, nodes)
        st <- c('A1', 'B1', 'C1', 'D1', 'D1', 'E1', 'F1')
        ed <- c('B1', 'C1', 'D1', 'F1', 'E1', 'G1', 'G1')
        for (i in 1:length(st)) {
            A[st[i], ed[i]] <- 1
        }
        dip_lay <- matrix(c(
                2.04, 3.15,
                2.97, 2.98,
                3.32, 2.70,
                4.00, 2.00,
                3.97, 1.50,
                5.45, 1.63,
                4.95, 1.21),
                nrow=7, byrow=TRUE)
        rownames(dip_lay) <- rownames(A)
        colnames(dip_lay) <- c("x", "y")
        class(A) <- "mcgraph"
        attr(A, "type") <- "constellation"
        attr(A, "mode") <- "directed"
        attr(A, "layout") <- dip_lay
    } else if (name == "virgin") {
        nodes <- 14
        A <- matrix(0, nrow=nodes, ncol=nodes)
        rownames(A) <- colnames(A) <- mcg.autonames(LETTERS, nodes)
        st <- c('A1', 'I1', 'G1', 'G1', 'G1', 'D1', 'H1', 'K1', 'N1', 'B1', 'F1', 'F1', 'L1', 'J1')
        ed <- c('I1', 'G1', 'D1', 'H1', 'F1', 'E1', 'K1', 'N1', 'B1', 'H1', 'L1', 'J1', 'C1', 'M1')
        for (i in 1:length(st)) {
            A[st[i], ed[i]] <- 1
        }
        vir_lay <- matrix(c(
                 2.22, 4.43,
                11.72, 5.46,
                -7.84, 5.34,
                 3.04, 5.63,
                 2.67, 6.23,
                 0.85, 5.16,
                 5.05, 5.24,
                 8.49, 5.38,
                 2.67, 4.84,
                -3.86, 4.91,
                10.52, 5.86,
                -1.99, 5.33,
                -7.49, 4.88,
                13.28, 5.76),
                nrow=14, byrow=TRUE)
        rownames(vir_lay) <- rownames(A)
        colnames(vir_lay) <- c("x", "y")
        class(A) <- "mcgraph"
        attr(A, "type") <- "constellation"
        attr(A, "mode") <- "directed"
        attr(A, "layout") <- vir_lay
    } else {
        stop("Error: Unknown constellation, known ones are: Lion, Dipper, Virgin!")
    }
    return(A)
}

#' @title Creation of graphs with a central node and outgoing bands
#'
#' @description `mcg.cross` creates a central node with several outgoind bands of same
#' length.
#'
#' @details This graph generation function will generated graphs with one centreal node and several chains of same length attached to the central node.
#' The final number of nodes will be bands x length + 1.
#' @param bands number of bands outgoing from the central node
#' @param length length of nodes number of nodes to create
#' @return graph of class mcgraph with adjacency matrix for an undirected graph
#' @examples
#' cross=mcg.cross(bands=8,length=5)
#' plot(cross)
#' @export

mcg.cross = function (bands=4,length=4) {
    nodes=bands*length+1
    mt=matrix(0,nrow=nodes,ncol=nodes)
    rownames(mt)=colnames(mt)=mcg.autonames(LETTERS,nodes)
    y=rep(1:length,bands)
    x=2:nodes
    start= -length
    for (i in 1:length(x)) {
        yr=y[i]
        if (y[i]==1) {
            start = start+length
        } else {
            yr=yr+start
        }
        mt[yr,x[i]]=1
        mt[x[i],yr]=1
    }
    class(mt)="mcgraph"
    attr(mt,"type")="cross"
    attr(mt,"mode")="undirected"
    return(mt)
}

#' @title Creation of circular graphs
#'
#' @description `mcg.circular` creates a circular graph with given number of nodes
#'
#' @details This function create a simple circular graph with a given number of nodes. All nodes in the graph have degree two.
#' @param nodes number of nodes to create in the graph
#' @return graph of class mcgraph with adjacency matrix for an undirected graph
#' @author Detlef Groth <dgroth@uni-potsdam.de>
#' @examples
#' circ=mcg.circular(12)
#' plot(circ,layout='circle')
#' @export

mcg.circular = function (nodes=8) {
    mt=matrix(0,nrow=nodes,ncol=nodes)
    rownames(mt)=colnames(mt)=mcg.autonames(LETTERS,nodes)
    for (row in 1:(nrow(mt)-1)) {
        mt[row,row+1]=1
        mt[row+1,row]=1
    }
    mt[nodes,1]=1
    mt[1,nodes]=1
    class(mt)="mcgraph"
    attr(mt,"type")="circular"
    attr(mt,"mode")="undirected"
    return(mt)
}
#' @title Creation of lattice like graphs
#' @description `mcg.lattice` creates a lattice graph of size `dim x dim` with some optional centralizing edges.
#' @details This function create a lattice graph of size `dim x dim`. Additionally some centralizing edges can be added connecting the angular outer nodes more to the center. The parameter `centralize = 1` will add four edges to all four corner nodes, `centralize = 2` will add eight edges.
#' @param dim dimensions of the lattice matrix, i.e. number of rows and columns
#' @param centralize Level of centralization, 0 no default no edges are added, 1 four edges are added to each corner edge and so on.
#' @return Graph of class mcgraph with adjacency matrix for an undirected graph.
#' @author Detlef Groth <dgroth@uni-potsdam.de>
#' @examples
#' lat=mcg.lattice(7,centralize=2)
#' plot(lat,layout='grid')
#' @export

mcg.lattice = function (dim=5,centralize=0) {
    mt=matrix(0,nrow=dim^2,ncol=dim^2)
    rownames(mt)=colnames(mt)=mcg.autonames(LETTERS,dim^2)
    for (i in 1:(dim-1)) {
        for (j in 0:(dim-1)) {
            mt[(i+(j*dim)+1),i+j*dim]=mt[i+(j*dim),(i+(j*dim)+1)]=1
        }
    }
    for (i in 1:dim) {
        for(j in 0:(dim-2)) {
            mt[(i+(j*dim)+dim),i+j*dim]=mt[i+(j*dim),(i+(j*dim)+dim)]=1
        }
    }
    if (centralize > 0 && centralize*2+1 > dim) {
        stop(paste("with centralize=",centralize,"dim must be at least ",centralize*2+1))

    }
    coords=c(1,2+dim,
             dim,2*dim-1,
             1+(dim*dim)-dim,dim*dim-(2*dim)+2,
             dim*dim,dim*dim-dim-1)
    if (centralize > 0) {
        for (i in 1:centralize) {
            mt[coords[1],coords[2]]=1
            mt[coords[3],coords[4]]=1
            mt[coords[5],coords[6]]=1
            mt[coords[7],coords[8]]=1
            ncoords=coords
            ncoords[1]=coords[2]
            ncoords[3]=coords[4]
            ncoords[5]=coords[6]
            ncoords[7]=coords[8]
            ncoords[2]=coords[2]+dim+1
            ncoords[4]=coords[4]+dim-1
            ncoords[6]=coords[6]-dim+1
            ncoords[8]=coords[8]-dim-1
            coords=ncoords
        }
    }
    mt[]  <- mt+t(mt)
    mt[mt>1]=1
    class(mt)="mcgraph"
    attr(mt,"type")="lattice"
    attr(mt,"mode")="undirected"
    return(mt)
}
#' @title Creation of random graphs
#'
#' @description `mcg.random` creates a random graph for given number of nodes and edges
#'
#' @details This function create a random graph whith a given number of nodes and edges.
#' All possible node pairs the same probability of being connected by an edge.
#' @param nodes Number of nodes in the graph.
#' @param edges Number of edges in the graph.
#' @return Graph of class mcgraph with adjacency matrix for an undirected graph.
#' @author Detlef Groth <dgroth@uni-potsdam.de>
#' @examples
#' rand=mcg.random(nodes=20,edges=30)
#' plot(rand,vertex.color=rainbow(20))
#' @export

mcg.random = function (nodes=26,edges=50) {
    mt=matrix(0,nrow=nodes,ncol=nodes)
    rownames(mt)=colnames(mt)=mcg.autonames(LETTERS,nodes)
    l=length(mt[upper.tri(mt)])
    zeros=rep(0,l)
    idx=sample(1:l,edges)
    zeros[idx]=1
    mt[upper.tri(mt)]=zeros
    mt[lower.tri(mt)]  <- t(mt)[lower.tri(mt)]
    class(mt)="mcgraph"
    attr(mt,"type")="random"
    attr(mt,"mode")="undirected"
    return(mt)
}

#' @title Creation of scale-free Barabasi-Albert graphs
#'
#' @description `mcg.barabasi` simplified mode to create a graph based on the Barabasi-Allbert algorithm
#'
#' @details This function create a random graph whith a given number of nodes.
#' Nodes are added subsequently connected by  defined number of edges to nodes added before in the graph.
#' @param nodes Number of nodes in the graph. Default: 26
#' @param m Number of edges to connect each added node to the graph. Default: 1
#' @return Graph of class mcgraph with adjacency matrix for an undirected graph.
#' @author Detlef Groth <dgroth@uni-potsdam.de>
#' @references Barabasi, A.-L. and Albert R. 1999. Emergence of scaling in random networks _Science_, 286 509-512
#' @references Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. http://igraph.org
#' @examples
#' bara=mcg.barabasi(nodes=20, m=1)
#' plot(bara,vertex.color=rainbow(20))
#' @export

mcg.barabasi = function (nodes=26, m=1) {
    # m number of edges to add at each node addition
    # second node can only have on edge
    mt=matrix(0,nrow=nodes,ncol=nodes)
    rownames(mt)=colnames(mt)=mcg.autonames(LETTERS, nodes)
    # MN Not working for higer number of nodes
    #mt['B1','A1']=1
    mt[2, 1]=1
    for (n in 3:ncol(mt)) {
        if (m==n) {
            sel=n-1
        } else {
            sel=m
        }
        idx=sample(1:(n-1),sel)
        mt[n,idx]=1
    }
    mt[upper.tri(mt)]  <- t(mt)[upper.tri(mt)]
    class(mt)="mcgraph"
    attr(mt,"type")="barabasi"
    attr(mt,"mode")="undirected"
    return(mt)
}

#' @title Creation of densely connected graph components
#'
#' @description `mcg.cluster` creates several graph components which are highly connected within each component.
#'
#' @details This function create a random graph whith a given number of nodes and a given cluster number.
#' Nodes are only connected to other nodes of the same cluster,
#' @param nodes Number of nodes in the graph. Default: 26
#' @param cluster Number of clusters, graph components. Default: 2
#' @inheritParams mcg.random
#' @return Graph of class mcgraph with adjacency matrix for an undirected graph.
#' @author Detlef Groth <dgroth@uni-potsdam.de>
#' @examples
#' clust=mcg.cluster(nodes=30,cluster=3)
#' plot(clust,vertex.color=rainbow(20))
#' @export
#'
mcg.cluster = function (nodes=26,cluster=2,edges=50) {
    mx=(nodes^2-nodes)/2
    if (mx < edges) {
        stop(paste("Error: number of edges exceeds number of possible edges which is",mx))
    }
    mt=matrix(0,nrow=nodes,ncol=nodes)
    rownames(mt)=colnames(mt)=mcg.autonames(LETTERS,nodes)
    sizes=rep(nodes %/% cluster,cluster)
    names(sizes)=1:cluster
    mod=nodes %% cluster ;
    if (mod > 0) {
        sizes=sizes+c(rep(1,mod),rep(0,cluster-mod))
    }
    edgs=rep(edges %/% cluster,cluster)
    names(edgs)=1:cluster
    mod=edges %% cluster ;
    if (mod > 0) {
        edgs=edgs+c(rep(1,mod),rep(0,cluster-mod))
    }
    cursor=1
    for (i in 1:cluster) {
        ang=as.matrix(mcg.angie(nodes=sizes[i],edges=edgs[i]))
        mt[cursor:(cursor+sizes[i]-1),cursor:(cursor+sizes[i]-1)]=ang
        cursor=cursor+sizes[i]
    }
    class(mt)="mcgraph"
    attr(mt,"type")="cluster"
    attr(mt,"mode")="undirected"
    return(mt)
}
#' @title Creation of small-world like random graphs within one component
#'
#' @description `mcg.angie` creates a single component graph with properties close to small world graphs.
#'
#' @details This function creates a random graph with a given number of nodes and a given number of edges by first building a tree like graph structure,
#' where nodes are stepwiese added to existing nodes in the graph. In the second step edges will be added which randomly connect two nodes in the graph.
#' @param nodes Number of nodes in the graph. Default: 26
#' @param edges Number of edges in the graph. Default: 40
#' @return Graph of class mcgraph with adjacency matrix for an undirected graph.
#' @author Detlef Groth <dgroth@uni-potsdam.de>
#' @references Fritz, Makejeva, Staub, Groth AA 2019.
#' @examples
#' angie=mcg.angie(nodes=26,edges=40)
#' plot(angie,vertex.color=rainbow(26))
#' @export

mcg.angie = function (nodes=26,edges=52) {
    mx=(nodes^2-nodes)/2
    if (mx < edges) {
        stop(paste("Error: number of edges exceeds number of possible edges which is",mx))
    }
    mt=matrix(0,nrow=nodes,ncol=nodes)
    rownames(mt)=colnames(mt)=mcg.autonames(LETTERS,nodes)
    mt[1,2]=1
    mt[2,1]=1
    for (e in 3:(edges+1)) {
        if (e <= nodes) {
            # Why not just colSums, which is way faster (~10x)? M.N.
            #idx=which(apply(mt,2,sum)>0)
            idx=which(colSums(mt) > 0)
            node=sample(idx,1)
            mt[e,node]=1
            mt[node,e]=1
        } else {
            while (TRUE) {
                nods=sample(1:nodes,2)
                if (mt[nods[1],nods[2]]==0) {
                    mt[nods[1],nods[2]]=1
                    mt[nods[2],nods[1]]=1
                    break
                }
            }

        }
    }
    class(mt)="mcgraph"
    attr(mt,"type")="angie"
    attr(mt,"mode")="undirected"
    return(mt)
}

#' @title Creation of graphs with a separate cluster of nodes with a central hub node
#'
#' @description `mcg.hubs` creates multiple graph components each having a central hub node.
#'
#' @details This function creates a graph with one to many components where every component has a central hub connected to all other nodes of the component.
#' @param nodes Number of nodes in the graph.
#' @param hubs Number of components each having one hub node.
#' @return Graph of class mcgraph with adjacency matrix for an undirected graph.
#' @author Detlef Groth <dgroth@uni-potsdam.de>
#' @examples
#' hubs=mcg.hubs(nodes=21,hubs=3)
#' plot(hubs,vertex.color=rainbow(21))
#' @export

mcg.hubs = function (nodes=8,hubs=2) {
    mt=matrix(0,nrow=nodes,ncol=nodes)
    rownames(mt)=colnames(mt)=mcg.autonames(LETTERS,nodes)
    mod=nodes %% hubs ;
    size = nodes %/% hubs;
    nhubs = (nodes-mod)/size ;
    cursor=0
    end=0
    start=1
    for (i in 1:nhubs) {
        if (mod > 0) {
            hsize = size+1;
            mod=mod-1
        } else {
            hsize = size
        }
        end=end+hsize
        for (j in (start+1):end) {
            mt[start,j]=1
            mt[j,start]=1
        }
        start=start+hsize
    }
    class(mt)="mcgraph"
    attr(mt,"type")="hubs"
    attr(mt,"mode")="undirected"
    return(mt)
}

#' @title Utility function to create name series
#'
#' @description `mcg.autonames` is a utility function to create sequences of names usable to nname nodes, edges, matrice rownames etc.
#'
#' @details This function simplifies the creation of names for a large set of nodes, edges, columns, rows etc. The user can give either a list of name prefixes or a single name prefix which will be used as prefix for the names.
#' The names are formatted with as much as required leading zeros.
#' @param nms name prefixes, can be one letter or a letter sequence.
#' @param n number of names to generate.
#' @return List of automatically generated names.
#' @author Detlef Groth <dgroth@uni-potsdam.de>
#' @examples
#' mcg.autonames("R",50)
#' mcg.autonames(LETTERS[1:3],20)
#' @export

mcg.autonames = function (nms,n) {
    # Arguments:
    #    nms: prefix or list of prefixes
    #    n:   number of names
    # Returns:
    #    vector of names of length n
    # Examples:
    #    autonames(c('A','B','C'), 10 )
    #  -> A1 B1 C1 A2 B2 C2 A3 B3 C3 A4
    #    autonames('R',12)
    #  -> R01 R02 .. R11 R12
    #
    # res=c()
    # if (length(nms)>1) {
    #     ln=n/length(nms)
    #     ni=floor(log10(ln)+1)
    #     ln=ln+1
    #     x=0
    #     for (i in 1:ln) {
    #         for (j in 1:(length(nms))) {
    #             nm=nms[j]
    #             format=paste(nm,"%0",ni,"i",sep="")
    #             res=c(res,sprintf(format,i))
    #         }
    #     }
    #     return(res[1:n])
    # } else {
    #     ni=floor(log10(n)+1)
    #     format=paste(nms,"%0",ni,"i",sep="")
    #     for (i in 1:n) {
    #         res=c(res,sprintf(format,i))
    #     }
    #     return(res[1:n])
    # }
    # Ommiting the for loops, leads to ~30x faster code. M.N.
    ln <- ceiling(n / length(nms))
    nms_tmp <- rep(nms, ln)[1:n]
    ptf <- rep(1:ln, each=length(nms))[1:n]
    frmt <- formatC(ptf, flag="0", width=log10(ln) + 1)
    return(paste0(nms_tmp, frmt))
}

#' @title Creation of special graphs
#'
#' @description `mcg.special` creates some very special directed graphs.
#'
#' @details This function just returns some interesting small graphs for playing and analysing.
#' @param gname Name of the graph to get. Possible names are currently: hub40, hub04, hub31, hub13, hub22,
#' hub30, hub03, hub21, hub12, nicolas, werner, wernerr, wernerc.
#' @return directed graph of class mcgraph.
#' @author Detlef Groth <dgroth@uni-potsdam.de>
#' @examples
#' g40=mcg.special("hub40");
#' plot(g40,vertex.color=rainbow(5))
#' nic=mcg.special("nicolas")
#' plot(nic,vertex.color=c("skyblue",rep("salmon", 4)),vertex.size=4,
#'      edge.width=3,arr.length=0.2,label.size=2,
#'      edge.color="black")
#' @export

mcg.special = function (gname) {
    if (gname == "nicolas") {
        mt=matrix(0,nrow=5,ncol=5)
        rownames(mt)=colnames(mt)=LETTERS[1:5]
        mt['A',c('B','E')]=1
        mt['B',c('C','E')]=1
        mt['C','D']=1
        mt['D',c('A','B')]=1
        mt['E','D']=1
        class(mt)="mcgraph"
        attr(mt,"type")="nicolas"
        attr(mt,"mode")="directed"
        lay=matrix(c(0,0, 0,1, 0.5,1.75, 1,1, 1,0),
                   ncol=2,byrow=TRUE)
        rownames(lay)=rownames(mt)
        colnames(lay)=c("x","y")
        attr(mt,'layout')=lay
        return(mt)
    } else if (gname == "werner") {
        mt=matrix(0,nrow=6,ncol=6)
        rownames(mt)=colnames(mt)=LETTERS[1:6]
        mt[c('A','B'),'C']=1
        mt['C','D']=1
        mt['D',c('E','F')]=1
        mt['E','F']=1
        class(mt)="mcgraph"
        attr(mt,"type")="werner"
        attr(mt,"mode")="directed"
        lay=matrix(c(1,2, 1,1, 2,1.5, 3,1.5, 4,2, 4,1),
                   ncol=2,byrow=TRUE)
        rownames(lay)=rownames(mt)
        colnames(lay)=c("x","y")
        attr(mt,'layout')=lay
        return(mt)
    } else if (gname == "wernerr") {
        mt=matrix(0,nrow=7,ncol=7)
        rownames(mt)=colnames(mt)=LETTERS[1:7]
        mt[c('A','B'),'C']=1
        mt['C','D']=1
        mt['D',c('E','F')]=1
        mt['E','G']=1
        mt['G','F']=1
        class(mt)="mcgraph"
        attr(mt,"type")="wernerr"
        attr(mt,"mode")="directed"
        lay=matrix(c(1,2, 1,1, 2,1.5, 3,1.5, 4,2, 4,1, 5,1.5), ncol=2,byrow=TRUE)
        rownames(lay)=rownames(mt)
        colnames(lay)=c("x","y")
        attr(mt,'layout')=lay
        return(mt)
    } else if (gname == "wernerc") {
        mt=matrix(0,nrow=8,ncol=8)
        rownames(mt)=colnames(mt)=LETTERS[1:8]
        mt[c('A','B'),'C']=1
        mt['C','D']=1
        mt['D',c('E','F')]=1
        mt['E','G']=1
        mt['G','H']=1
        mt['H','F']=1
        class(mt)="mcgraph"
        attr(mt,"type")="wernerc"
        attr(mt,"mode")="directed"
        lay=matrix(c(1,2, 1,1, 2,1.5, 3,1.5, 4,2, 4,1, 5,2, 5,1), ncol=2,byrow=TRUE)
        rownames(lay)=rownames(mt)
        colnames(lay)=c("x","y")
        attr(mt,'layout')=lay
        return(mt)
    } else if (gname == "mini") {
        mt=matrix(0,nrow=4,ncol=4)
        rownames(mt)=colnames(mt)=LETTERS[1:4]
        mt['A','B']=1;
        mt['B','C']=1;
        mt['D','C']=1;
        class(mt)="mcgraph"
        attr(mt,"type")="mini"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "hub40") {
        mt=matrix(0,nrow=5,ncol=5)
        rownames(mt)=colnames(mt)=LETTERS[1:5]
        mt['A',c('B','C','D','E')]=1
        class(mt)="mcgraph"
        attr(mt,"type")="hub40"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "hub04") {
        mt=matrix(0,nrow=5,ncol=5)
        rownames(mt)=colnames(mt)=LETTERS[1:5]
        mt[c('B','C','D','E'),'A']=1
        class(mt)="mcgraph"
        attr(mt,"type")="hub04"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "hub31") {
        mt=matrix(0,nrow=5,ncol=5)
        rownames(mt)=colnames(mt)=LETTERS[1:5]
        mt['A',c('B','C','D')]=1
        mt['E','A']=1
        class(mt)="mcgraph"
        attr(mt,"type")="hub31"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "hub13") {
        mt=matrix(0,nrow=5,ncol=5)
        rownames(mt)=colnames(mt)=LETTERS[1:5]
        mt[c('B','C','D'),'A']=1
        mt['A','E']=1
        class(mt)="mcgraph"
        attr(mt,"type")="hub13"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "hub22") {
        mt=matrix(0,nrow=5,ncol=5)
        rownames(mt)=colnames(mt)=LETTERS[1:5]
        mt[c('B','C'),'A']=1
        mt['A',c('D','E')]=1
        class(mt)="mcgraph"
        attr(mt,"type")="hub22"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "hub30") {
        mt=matrix(0,nrow=4,ncol=4)
        rownames(mt)=colnames(mt)=LETTERS[1:4]
        mt['A',c('B','C','D')]=1
        class(mt)="mcgraph"
        attr(mt,"type")="hub30"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "hub03") {
        mt=matrix(0,nrow=4,ncol=4)
        rownames(mt)=colnames(mt)=LETTERS[1:4]
        mt[c('B','C','D'),'A']=1
        class(mt)="mcgraph"
        attr(mt,"type")="hub03"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "hub21") {
        mt=matrix(0,nrow=4,ncol=4)
        rownames(mt)=colnames(mt)=LETTERS[1:4]
        mt['A',c('B','C')]=1
        mt['D','A']=1
        class(mt)="mcgraph"
        attr(mt,"type")="hub21A"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "hub21B") {
        mt=matrix(0,nrow=4,ncol=4)
        rownames(mt)=colnames(mt)=LETTERS[1:4]
        mt['A',c('B','C')]=1
        mt['D','A']=1
        mt['C','B']=1
        class(mt)="mcgraph"
        attr(mt,"type")="hub21B"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "hub12") {
        mt=matrix(0,nrow=4,ncol=4)
        rownames(mt)=colnames(mt)=LETTERS[1:4]
        mt[c('B','C'),'A']=1
        mt['A','D']=1
        class(mt)="mcgraph"
        attr(mt,"type")="hub12"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "3A") {
        mt=matrix(0,nrow=3,ncol=3)
        rownames(mt)=colnames(mt)=LETTERS[1:3]
        mt['A',c('B','C')]=1
        class(mt)="mcgraph"
        attr(mt,"type")="3A"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "3B") {
        mt=matrix(0,nrow=3,ncol=3)
        rownames(mt)=colnames(mt)=LETTERS[1:3]
        mt[c('B','C'),'A']=1
        mt['B','C']=1
        class(mt)="mcgraph"
        attr(mt,"type")="3B"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "3C") {
        mt=matrix(0,nrow=4,ncol=4)
        rownames(mt)=colnames(mt)=LETTERS[1:4]
        mt['A','B']=1
        mt['B','C']=1
        mt['B','D']=1
        class(mt)="mcgraph"
        attr(mt,"type")="3C"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "triA") {
        mt=matrix(0,nrow=4,ncol=4)
        rownames(mt)=colnames(mt)=LETTERS[1:4]
        mt['A','B']=1
        mt['B','C']=1
        class(mt)="mcgraph"
        attr(mt,"type")="triA"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "triB") {
        mt=matrix(0,nrow=4,ncol=4)
        rownames(mt)=colnames(mt)=LETTERS[1:4]
        mt['C','B']=1
        mt['B','A']=1
        class(mt)="mcgraph"
        attr(mt,"type")="triB"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "triC") {
        mt=matrix(0,nrow=4,ncol=4)
        rownames(mt)=colnames(mt)=LETTERS[1:4]
        mt['A','B']=1
        mt['B','C']=1
        mt['A','D']=1
        mt['D','C']=1
        class(mt)="mcgraph"
        attr(mt,"type")="triC"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "triD") {
        mt=matrix(0,nrow=4,ncol=4)
        rownames(mt)=colnames(mt)=LETTERS[1:4]
        mt['A','B']=1
        mt['B','C']=1
        mt['A','C']=1
        class(mt)="mcgraph"
        attr(mt,"type")="triD"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "triE") {
        mt=matrix(0,nrow=4,ncol=4)
        rownames(mt)=colnames(mt)=LETTERS[1:4]
        mt['A','B']=1
        mt['B','C']=1
        mt['D','C']=1
        class(mt)="mcgraph"
        attr(mt,"type")="triE"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "triF") {
        mt=matrix(0,nrow=4,ncol=4)
        rownames(mt)=colnames(mt)=LETTERS[1:4]
        mt['A','B']=1
        mt['B','C']=1
        mt['D','A']=1
        class(mt)="mcgraph"
        attr(mt,"type")="triF"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "triG") {
        mt=matrix(0,nrow=4,ncol=4)
        rownames(mt)=colnames(mt)=LETTERS[1:4]
        mt['A','B']=1
        mt['B','C']=1
        mt['D','B']=1
        class(mt)="mcgraph"
        attr(mt,"type")="triG"
        attr(mt,"mode")="directed"
        return(mt)
    } else if (gname == "tricky") {
        mt=matrix(0,nrow=9,ncol=9)
        rownames(mt)=colnames(mt)=LETTERS[1:9]
        mt['A','B']=1
        mt['A','C']=1
        mt['B','D']=1
        mt['C','I']=1
        mt['D','E']=1
        mt['E','F']=1
        mt['F','G']=1
        mt['G','H']=1
        mt['H','I']=1
        mt['I','D']=1
        class(mt)="mcgraph"
        attr(mt,"type")="tricky"
        attr(mt,"mode")="directed"
        lay=matrix(c(1,2.5, 2,3, 2,2, 3,3, 4,3.5, 5,3,
                     5,2, 4,1.5, 3,2), ncol=2,byrow=TRUE)
        rownames(lay)=rownames(mt)
        colnames(lay)=c("x","y")
        attr(mt,'layout')=lay
        return(mt)
    } else if (gname == "caputh") {
        mt=matrix(0,nrow=30,ncol=30)
        rownames(mt)=colnames(mt)=c(LETTERS[1:14],
                                    paste(LETTERS[1:4],1,sep=""),
                                    paste(LETTERS[1:4],2,sep=""),
                                    paste(LETTERS[1:4],3,sep=""),
                                    paste(LETTERS[1:4],4,sep=""))
        mt['A','B']=1;
        mt['B','C']=1;
        mt['C','D']=1
        mt['D','E']=1;
        mt['E','F']=1;
        mt['F','G']=1;
        mt['G','H']=1
        mt['H','I']=1;
        mt['I','J']=1;
        mt['J','K']=1
        mt['K','L']=1
        mt['L','M']=1
        mt['M','N']=1
        mt['C','A1']=1
        mt['A1','B1']=1
        mt['B1','C1']=1
        mt['D1','A1']=1
        mt['F','A2']=1
        mt['A2','B2']=1
        mt['B2','C2']=1
        mt['D2','B2']=1
        mt['I','A3']=1
        mt['A3','B3']=1
        mt['B3','C3']=1
        mt['D3','C3']=1
        mt['L','A4']=1
        mt['A4','B4']=1
        mt['B4','C4']=1
        mt['A4','D4']=1
        mt['D4','C4']=1
        class(mt)="mcgraph"
        attr(mt,"type")="caputh"
        attr(mt,"mode")="directed"
        attr(mt,'layout')=matrix(c(rep(1,14),2,3,4,3,2,3,4,3,2,3,4,3,2,3,4,3,
            15:2, 13,14,13,12, 10,11,10,9, 7,8,7,6, 4,5,4,3),ncol=2,byrow=FALSE)
        return(mt)
    } else {
       stop("unkown graph name")
    }
}
#' @title mcg.corrplot simple display of correlation matrices
#' @description `mcg.corrplot` provides a simple display of correlation matrices using filled circles to represent the strength of correlation.
#' @param mt correlation matrix or any other numeric matrix as input. Values must be in the range of -1, 1.
#' @param text.lower should the coefficients be plottet in the lower plot triangle. Default: FALSE
#' @param text.upper should the coefficients be plottet in the upper plot triangle. Default: FALSE
#' @param pch.plus plotting character for positive correlations. Default: 19
#' @param pch.minus plotting character for negative correlations, a good possibility is 18. Default: 19
#' @param xtext text added at the bottom of the plot below of the last value line.
#' @param cex relative size of the plotting symbols on the plot. Default: 1
#' @param cex.coeff enlargement coefficient for the r coefficients. Default: 1
#' @param cex.names enlargment value for row and column names. Default: 1
#' @param \dots Any other argument will be passed to the generic plot function.
#' @return nothing.
#' @details This is a small plot utility to display correlation matrices and other matrices with numeric values in the range of -1, 1.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  ang=mcg.angie(nodes=12,edges=15)
#'  anu=mcg.u2d(ang)
#'  and=mcg.graph2data(anu)
#'  mcg.corrplot(cor(t(and)),text.lower=TRUE)
#'  pca=prcomp(scale(t(and))
#'  mcg.corrplot(pca$rotation)
#'  }
#' }
#' @rdname mcg.corrplot
#' @export

mcg.corrplot = function (mt,text.lower=FALSE,text.upper=FALSE,pch.minus=19,
                         pch.plus=19,xtext=NULL,cex=1,cex.names=1,
                         cex.coeff=1,...) {
    plot(1,type="n",xlab="",ylab="",axes=FALSE,
         xlim=c(-0.5,ncol(mt)+0.5),ylim=c(nrow(mt)+0.5,0),...)

    #cex=cex*0.9*(10/nrow(mt))
    cex.names=cex.names*0.9*(10/nrow(mt))
    #if (cex>0.9) {
    #    cex=1
    #}
    # change
    text(1:ncol(mt),0.25,colnames(mt),cex=cex.names)
    if (length(xtext)>0) {
        text(1:ncol(mt),nrow(mt)+0.75,xtext,cex=cex.names)
    }
    # change
    text(0,1:nrow(mt),rownames(mt),cex=cex.names)
    cols=paste("#DD3333",rev(c(15,30, 45, 60, 75, 90, "AA","BB","CC","DD")),sep="")
    cols=c(cols,paste("#3333DD",c(15,30, 45, 60, 75, 90, "AA","BB","CC","DD"),sep=""))
    breaks=seq(-1,1,by=0.1)
    sym=identical(rownames(mt),colnames(mt))
    cex=3*(10/nrow(mt))*cex
    #if (cex>5) {
    #    cex=5
    #}
    for (i in 1:nrow(mt)) {
        for (j in 1:ncol(mt)) {
            if (sym & i == j) {
                next
            }
            #cex=abs(mt[i,j])*2
            coli=cut(mt[i,j],breaks=breaks,labels=1:20)
            pch=19
            if (is.na(mt[i,j])) {
                pch=NA
            } else if (mt[i,j]< 0) {
                pch=pch.minus
            } else {
                pch=pch.plus
            }
            if (pch==17) {
                points(i,j,pch=pch,cex=cex*0.8,col=cols[coli])
            } else if (!is.na(pch)) {
                points(i,j,pch=pch,cex=cex,col=cols[coli])
            }
            if (i == j & !sym & text.lower) {
                points(i,j,pch=19,cex=cex*1.1,col='white')
                text(i,j,sprintf("%.2f",mt[i,j]),cex=cex.coeff)
            } else if (i < j & text.lower) {
                points(i,j,pch=19,cex=cex*1.1,col='white')
                text(i,j,sprintf("%.2f",mt[i,j]),cex=cex.coeff)
            } else if (i > j & text.upper) {
                points(i,j,pch=19,cex=cex*1.1,col='white')
                text(i,j,sprintf("%.2f",mt[i,j]),cex=cex.coeff)
            }

        }
    }
}
#'
#' @title Evaluate the time required for a given R expression
#' @description The function can be used to measure the required computation time for a given R expression
#'   which is executed one or more times.
#' @param expr a valid R expression
#' @param n how often should the expression be executed, default: 5
#' @return a vector with n values in seconds the time required to evaluate the expression
#' @examples
#' data(iris)
#' mcg.timeit(expression(hclust(dist(iris[,1:4]))))
#' data(swiss)
#' mcg.timeit(expression(mcg.lvs(swiss, rs=0.1, output='mcgraph')))
#' mean(mcg.timeit(expression(mcg.lvs(swiss, rs=0.1, output='mcgraph'))))
#' mean(mcg.timeit(expression(mcg.lvs(swiss, rs=0.1, output='mcgraph',code="R"))))
#' @rdname mcg.timeit
#' @export
#' @author Detlef Groth <email: dgroth@uni-potsdam.de>
#'

mcg.timeit = function (expr,n=5) {
    if (!is.expression(expr)) {
        stop("Error: given expression must be created with expression function such as: \nexpression(dist(iris[,1:4])) !")
    }
    res=c()
    for (i in 1:n) {
        t1=Sys.time();
        eval(expr);
        res=c(res,as.numeric(difftime(Sys.time(),t1,units="sec")))
    }
    return(res)
}

