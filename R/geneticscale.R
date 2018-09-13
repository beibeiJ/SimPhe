##' Genetic scales of epistatic model (Cockerham model) based on F2 population.
##'
##' Calculate the genetic scale for a given allele combination of two SNPs.
##' There are 9 genotypes in an F2 population, so we need 8 genetic parameters
##' to give a complete description of the values for the 9 genotypes. Under the assumption
##' of Hardy-Weinberg and linkage equilibrium, Cockerham (1954)'s orthogonal
##' partition of genetic variance leads to the definition of the genotypic
##' value \eqn{G_{ij}}
##'
##' \eqn{G_{ij}=\beta_{0}+\sum_{t=1}^{8}\beta_{G_{w_{t}}}w_{tij}}
##'
##' by eight orthogonal scales or contrasts \eqn{w_{t}}'s, t in (1, 2, \eqn{\cdots}{...}, 8). Four are marginal
##' scales and four are interaction scales. Marginal scales (defined by Model I
##' for an F2 population) are called linear and quadratic scales (additive and
##' dominance scales in genetic terms). Correspondingly, the interaction scales
##' are
##' \describe{
##'   \item{w1}{additive for locus A;}
##'   \item{w2}{dominance for locus A;}
##'   \item{w3}{additive for locus B;}
##'   \item{w4}{dominance for locus B;}
##'   \item{w5 (= w1 \eqn{\times}{*} w3)}{linear    \eqn{\times}{*} linear, additive  \eqn{\times}{*} additive;}
##'   \item{w6 (= w1 \eqn{\times}{*} w4)}{linear    \eqn{\times}{*} quadratic, additive  \eqn{\times}{*} dominance;}
##'   \item{w7 (= w2 \eqn{\times}{*} w3)}{quadratic \eqn{\times}{*} linear, dominance \eqn{\times}{*} additive;}
##'   \item{w8 (= w2 \eqn{\times}{*} w4)}{quadratic \eqn{\times}{*} quadratic, dominance \eqn{\times}{*} dominance.}
##' }
##' SNPs are encoded by (0, 1, 2):
##' \describe{
##'   \item{0}{means homozygous with major alleles;}
##'   \item{1}{means heterozygote;}
##'   \item{2}{means homozygous with minor alleles}
##' }
##' E.g., the SNPs are encoded as 0: AA, 1: AG, and 2: GG,  where 'A' represents the major allele and 'G' the minor allele.
##' @title Genetic scales of epistatic model
##' @export
##' @param SNPA encoded alleles for first SNP.
##' @param SNPB encoded alleles for second SNP.
##' @return a vector of genetic scales.
##' @author Benno Pütz \email{puetz@@mpipsykl.mpg.de} and Beibei Jiang \email{beibei_jiang@@psych.mpg.de}
##' @examples
##' genetic.scale(SNPA = 1, SNPB = 0)
genetic.scale <- function(SNPA = 0,
                          SNPB = 0){
    fixSNP <- function(snp) min(2, max(0, trunc(snp+0.5))) # don't want IEEE rounding of round()
    SNPA <- fixSNP(SNPA)
    SNPB <- fixSNP(SNPB)
## message(c(SNPA, SNPB))
    w1 <- 1 - SNPA
    w2 <- ifelse(SNPA==1, 0.5, -0.5)
    w3 <- 1 - SNPB
    w4 <- ifelse(SNPB==1, 0.5, -0.5)
    return(c(w1, w2, w3, w4,
             w1 * w3,
             w1 * w4,
             w2 * w3,
             w2 * w4))
}





