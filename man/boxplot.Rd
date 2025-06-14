% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boxplot.r
\name{boxplot}
\alias{boxplot}
\title{Boxplot}
\usage{
boxplot(
  dat,
  factor,
  outline = FALSE,
  title = "Boxplot",
  legendtitle = "TREATMENT",
  outlier = T,
  coefiqr = 1.5,
  ggplot = FALSE
)
}
\arguments{
\item{dat}{matrix numeric}

\item{factor}{factor}

\item{outline}{logical display outliers FALSE by default}

\item{title}{character}

\item{legendtitle}{character}

\item{outlier}{boolean}

\item{coefiqr}{numeric}

\item{ggplot}{logical use graphics library or ggplot FALSE by default}
}
\value{
plot
}
\description{
Boxplot
}
\details{
To make boxplot from matrix.
}
\examples{
# not run
# mat1 \%>\% boxplot( factor = sif1$F3 )
}
\author{
Florent Dumont \href{mailto:florent.dumont@universite-paris-saclay.fr}{florent.dumont@universite-paris-saclay.fr}
}
