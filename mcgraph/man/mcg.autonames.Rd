% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcgraph.R
\name{mcg.autonames}
\alias{mcg.autonames}
\title{Utility function to create name series}
\usage{
mcg.autonames(nms, n)
}
\arguments{
\item{nms}{name prefixes, can be one letter or a letter sequence.}

\item{n}{number of names to generate.}
}
\value{
List of automatically generated names.
}
\description{
`mcg.autonames` is a utility function to create sequences of names usable to nname nodes, edges, matrice rownames etc.
}
\details{
This function simplifies the creation of names for a large set of nodes, edges, columns, rows etc. The user can give either a list of name prefixes or a single name prefix which will be used as prefix for the names.
The names are formatted with as much as required leading zeros.
}
\examples{
mcg.autonames("R",50)
mcg.autonames(LETTERS[1:3],20)
}
\author{
Detlef Groth <dgroth@uni-potsdam.de>
}
