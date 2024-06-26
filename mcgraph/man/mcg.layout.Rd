% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcgraph.R
\name{mcg.layout}
\alias{mcg.layout}
\title{Layout a graph for plotting}
\usage{
mcg.layout(A, mode = "sam", noise = FALSE, ...)
}
\arguments{
\item{A}{Adjacency matrix}

\item{mode}{Character string for the plotting algorithm. Available are 'mds', 'sam', 'frg', 'circle' and 'grid', Default: sam}

\item{noise}{Boolean. Should noise be added to the coordinates. Useful if vertices overlap.}

\item{...}{Other parameters concerning individual layout functions. See 'mcg.layout.frg'.}
}
\value{
matrix with two columns for x and y coordinates.
}
\description{
Layout a graph for plotting
}
\examples{
ang=mcg.angie(nodes=12,edges=16)
xy=mcg.layout(ang,mode='mds')
plot(ang,layout=xy)
}
\author{
Detlef Groth <email: dgroth@uni-potsdam.de>
}
