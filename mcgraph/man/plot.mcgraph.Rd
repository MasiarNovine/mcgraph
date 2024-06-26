% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcgraph.R
\name{plot.mcgraph}
\alias{plot.mcgraph}
\title{Visualizing graph objects}
\usage{
\method{plot}{mcgraph}(
  x,
  layout = NULL,
  noise = FALSE,
  vertex.size = 1,
  vertex.length = 0.5,
  vertex.color = "orange",
  vertex.in.color = "light blue",
  vertex.out.color = "salmon",
  vertex.border.color = "black",
  vertex.symbol = "circle",
  edge.width = 1.5,
  edge.color = "black",
  edge.neg.color = "red",
  edge.lty = "solid",
  arr.length = 0.08,
  arr.angle = 12,
  arr.mode = "aligned",
  show.labels = TRUE,
  label.size = 1,
  label.font = NULL,
  label.color = NULL,
  label.adj = NULL,
  label.pos = NULL,
  label.offset = 0.5,
  label.text = NULL,
  label.vfont = NULL,
  bw = FALSE,
  interactive = FALSE,
  ...
)
}
\arguments{
\item{x}{Object of \code{mcgraph} or matrix class}

\item{layout}{Either a character string defining how to plot the graph with
possible values are \code{"mds"}, \code{"sam"}, \code{"frg"}, \code{"circle"}
and \code{"grid"} or a matrix with x and y coordinates for each node. Default: \code{"frg"}}

\item{noise}{Should be noise added to the layout? Sometimes useful if nodes are too close. Default: \code{FALSE}}

\item{vertex.size}{Size of nodes. If \code{vertex.symbol="rectangle"}, width of rectangles.  Default: \code{1}}

\item{vertex.length}{Height of nodes, if \code{vertex.symbol="rectangle"}. Default: \code{0.5}}

\item{vertex.color}{Color of nodes. Default: \code{"grey80"}}

\item{vertex.in.color}{Color of in-nodes in directed graphs. Default: \code{"light blue"}}

\item{vertex.out.color}{Color of out-nodes in directed graphs. Default: \code{"salmon"}}

\item{vertex.border.color}{Color of the border of nodes. Default: \code{"black"}}

\item{vertex.symbol}{Shape of node symbol. Possible values: \code{"circle"}, \code{"rectangle"}. Default: \code{"circle"}}

\item{edge.width}{Thickness of edges. Default: \code{1.5}}

\item{edge.color}{Color of positive edges. Default: "\code{grey50}"}

\item{edge.neg.color}{Color of negative edges. Default: \code{"red"}}

\item{edge.lty}{Linetype of edges. Default: \code{"solid"}}

\item{arr.length}{Length of arrows. Default: \code{0.08}}

\item{arr.angle}{Angle of arrows. Default: \code{12}}

\item{arr.mode}{Specifies where the arrow heads are plotted. Either \code{"midpoint"} or \code{"aligned"}. Default: \code{"aligned"}}

\item{show.labels}{Should node labels be printed? Default: \code{TRUE}}

\item{label.size}{Size of the node text labels. Default: \code{1}}

\item{label.font}{Font type to be used. \code{1} = normal, \code{2} = bold, \code{3} = italic, \code{4} = bold-italic. Default: \code{NULL}}

\item{label.color}{Color of label. Default: \code{NULL}}

\item{label.adj}{Justification of labels. \code{O} for left/bottom, \code{1} for rigth/top, \code{0.5} for centered. Default: \code{NULL}}

\item{label.pos}{Position specifier. \code{1} = below, \code{2} = left, \code{3} = above, \code{1} = right. Default: \code{NULL}}

\item{label.offset}{Distance of label from specified coordinate if \code{label.pos} is specified. Default: \code{0.5}}

\item{label.text}{Use costum labels instead of automatic ones. Default: \code{NULL}}

\item{label.vfont}{Null for default font family or character vector of length 2 for Hershey vector fonts. Default: \code{NULL}}

\item{bw}{Boolean Should only black and white color schemes used? Default: \code{FALSE}}

\item{interactive}{Switch into interactive mode where you can click in the graph
and move nodes with two clicks, first selecting the node, second
click gives thehe new coordinates for the node. Default: \code{FALSE}}

\item{...}{Other graphical parameters.}
}
\value{
Returns the layout of the plotted network (invisible).
}
\description{
Plot function to visualize \code{mcgraph} objects.
}
\details{
This is the basic plot functionality of \code{mcgraph} for plotting
         undirected and directed graphs. It is possible to give a layout
         for the graph by specifying the \code{layout} argument by a string
         or a two column matrix with the \code{x} and \code{y} coordinates
         for each node. The string must be one of the possible layout
         types given in \code{mcg.layout}, i.e. \code{mds} (multi-dimensional scaling), \code{sam}, \code{frg}
         (Fruchterman-Reingold force-directed layout), \code{circle} and \code{grid} (see \code{?mcg.layout} for more details).
         By default, the \code{frg} layout is used.

         If no \code{dimnames} are given, the nodes will be named by \code{mcg.autonames}.
         It is also possible to omit node labels completely by setting the \code{show.labels}
         to \code{FALSE} or to give costum \code{labels} in form of a string vector.

         Next to specifying basic properties like the color and the size for nodes, edges and the labels,
         there is also an \code{interactive} mode for modifying the layout by hand and saving the layout.
         Label properties can be specified via corresponding arguments used in the \code{text} function of the traditional \code{graphics}.

         The function offers basic plotting functionalities. For more sophisticated use-cases, the \code{igraph}
         package by G. Csardi might be more appropriate.
}
\examples{
band=mcg.band(16)
mcb=mcg.graph(band,n=200,input='A1')
plot(mcb,layout='circle')

# Larger network
mca <- mcg.angie(200, 320)
plot(mcb,layout='circle')
plot(mcb,layout='circle')
}
\author{
: Detlef Groth <email: dgroth@uni-potsdam.de>, Masiar Novine <email: masiar.novine@gmail.com>
}
