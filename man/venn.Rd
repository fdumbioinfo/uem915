% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/venn.r
\name{venn}
\alias{venn}
\title{Venn diagramm}
\usage{
venn(
  list = NULL,
  listnames = NULL,
  returnlist = F,
  title = "Venn Diagram",
  plot = T,
  export = F,
  path = ".",
  dirname = "venn"
)
}
\arguments{
\item{list}{list of 2 , 3 or 4 character vector or list of two data.frame to compare}

\item{listnames}{character list names to display on graph}

\item{returnlist}{logical}

\item{title}{character title to display on graph}

\item{plot}{logical to display the plot or not}

\item{export}{logical export list in file}

\item{path}{character}

\item{dirname}{character name of the directory created when export = T}
}
\value{
venn plot and new lists generated by venn.
}
\description{
To make a Venn diagramm of 2, 3 or 4 lists
}
\details{
until 4 list
}
\examples{
# not run
# library(magrittr)
# list(
#   c(letters[6:20] , letters[25] ) ,
#   letters[1:15] ,
#   c( letters[2:5] , letters[8:23] ) ) \%>\% venn
}
\author{
Florent Dumont \href{mailto:florent.dumont@universite-paris-saclay.fr}{florent.dumont@universite-paris-saclay.fr}
}
