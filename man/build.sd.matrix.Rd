% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlatedpheno.R
\name{build.sd.matrix}
\alias{build.sd.matrix}
\title{build a matrix with standard deviation}
\usage{
build.sd.matrix(x, margin = 2, ...)
}
\arguments{
\item{x}{a matrix or dataframe.}

\item{margin}{an integer specifying the dimension that standard
deviation will be computed for.
1 indicates rows, 2 indicates columns. Default is 2.}

\item{...}{not used}
}
\value{
a matrix with the standard deviation.
    If the number of dimensions is 2, then it is
    \eqn{\left(\begin{array}{cc} \sigma_{1} & 0 \\ 0 & \sigma_{2} \end{array}\right)}.
}
\description{
Build a matrix with standard deviation.
}
\examples{
x1 <- rnorm(4000, mean =  5, sd = 10)
x2 <- rnorm(4000, mean = 10, sd = 30)
x <- matrix(cbind(x1, x2), ncol = 2)
build.sd.matrix(x)
}
\author{
Beibei Jiang \email{beibei_jiang@psych.mpg.de} and
        Benno Pütz \email{puetz@psych.mpg.de}
}
