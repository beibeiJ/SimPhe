% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readpars.R
\name{read.simu.pars}
\alias{read.simu.pars}
\title{Read parameters}
\usage{
read.simu.pars(file = NULL, ...)
}
\arguments{
\item{file}{the file of parameters settings}

\item{...}{not used}
}
\value{
a list of simulation parameters. One element per block of \code{file}.
Each element is a dataframe with SNP names and model parameters
}
\description{
Read file specifying the simulation parameters.
}
\details{
File format: please follow the example of the \code{simupars.txt} file
found in the inst/extdata/ directory of the package (run \cr
\code{system.file("extdata", "simupars.txt", package="SimPhe")} \cr
to see the full path),
blank lines are ignored.
The file consists of three or four blocks for each phenotype:
mean, main, epistasis, and (optionally) heritability.
Each block is started by a line of the form '[blockname]' followed by the
parameters for the block, e.g., for thefirst phenotype,
\describe{
  \item{[P1mean]}{\describe{
                \item{mean}{\eqn{\beta_{0}}{\beta_0}: coefficient parameter of "basic" genetic effects in
                           \eqn{G_{ij}=\beta_{0}+\sum_{t=1}^{8}\beta_{G_{w_{t}}}w_{tij}}{Gij = \beta_0 + \sum(\beta_Gwt*wtij), t in (1, 2, ..., 8)}. }
                }}
  \item{[P1main]}{\describe{
                \item{SNP}{SNP name}
                \item{additive}{coefficient of additive effect}
                \item{dominance}{coefficient of dominance effect}
                }}
  \item{[P1epistasis]}{\describe{
                    \item{SNPA}{first SNP name}
                    \item{SNPB}{second SNP name}
                    \item{additive_additive}{coefficient for additive-additive interaction}
                    \item{additive_dominance}{coefficient for additive-dominance interaction}
                    \item{dominance_additive}{coefficient for dominance-additive interaction}
                    \item{dominance_dominance}{coefficient for dominance-dominance interaction}
                    }}
  \item{[P1heritability]}{\describe{
                        \item{heritability}{expected heritability}
                        }}
}
For each block the expected columns and the respective menaings are given.
Similar blocks need to be provided for the other phenotype(s): "[P2mean]",
"[P2main]", "[P2epistasis]", and so on.
}
\examples{
# simulation parameters:
fpar.path <- system.file("extdata", "simupars.txt", package="SimPhe")


# pars <- read.simu.pars(fpar.path)
}
\author{
Beibei Jiang \email{beibei_jiang@psych.mpg.de} and
        Benno Pütz \email{puetz@psych.mpg.de}
}
