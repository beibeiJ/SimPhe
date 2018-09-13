##' Calculate the gene(allele) frequency for each of the SNPs.
##' @title Calculate the gene(allele) frequency
##' @export
##' @param geno a dataframe of genotype data: columns are the SNPs; lines are indviduals.
##' @param epi.pars a data.frame or a matrix containing the parameter information
##'  for epistatic effect:
##' additive  \eqn{\times}{*} additive,
##' additive  \eqn{\times}{*} dominance,
##' dominance \eqn{\times}{*} additive, and
##' dominance \eqn{\times}{*} dominance.
##' @param ... not used
##' @return a dataframe with allele frequencies (major and minor).
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de}
##' @examples
##' # genotype file: rows are individuals and columns are SNPs
##' fgeno.path <- system.file("extdata", "10SNP.txt", package="SimPhe")
##'
##' # get genotype
##' geno <- read.geno(fgeno.path, ftype = "snp.head")
##'
##' get.freq(geno, epistasis.pars)
get.freq <- function(geno,
                     epi.pars,
                     ...){
  snps <- c(epi.pars[, "SNPA"], epi.pars[, "SNPB"])
  snps <- snps[!duplicated(snps)] # delete duplicates
  freq <- data.frame(matrix(nrow = length(snps), ncol = 3), stringsAsFactors = FALSE)
  freq[, 1] <- snps
  for(i in 1:length(snps))
   freq[i, -1] <- matrix(count.allele(geno[, snps[i]]), nrow = 1)
  colnames(freq) <- c("SNP", "major.frequency", "minor.frequency")
  return(freq)
}


##' Count major and minor allele frequencies.
##' @title Count allele frequencies
##' @export
##' @param x a vector of single SNP information (minor allele count for genotype).
##' @param ... not used.
##' @return a vector with major and minor allele frequency.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de} and
##'         Benno Pütz \email{puetz@@psych.mpg.de}
##' @examples
##' maf <- 0.1
##' x <- sample(0:2, 1000, replace = TRUE, prob = c((1-maf)^2, 2*(1-maf)*maf, maf^2))
##' table(x)
##' count.allele(x)
count.allele <- function(x,
                         ...){
  ncls <- c(table(x), rep(0))[1:3]
  N <- length(x) * 2
  # maj.freq <- (ncls["0"] * 2 + ncls["1"]) / (length(x) * 2)
  # min.freq <- (ncls["2"] * 2 + ncls["1"]) / (length(x) * 2)
  # return(c(maj.freq, min.freq))
  w <- 0:2
  major <- ncls %*% rev(w)
  minor <- ncls %*% w
  return(c(major,minor)/N)
 }
