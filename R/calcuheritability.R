##' Calculate heritability
##' @title Calculate heritability
##' @export
##' @param gene.coef a vector of 8 genetic parameters.
##'                  Each element includes 8 genetic parameters (regression
##'                  coefficient)
##'                  \eqn{\beta_{Gw_{t}}}, t in (1, 2, \eqn{\cdots}{...}, 8).
##' @param freq a dataframe with the allele frequencies.
##' @param noise.var variance of noise to generate the random noise.
##' @param Dskim the coefficient of linkage disequilibrium. Default is 0 (no LD).
##' @param ... not used.
##' @return heritability.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de}

calc.herit <- function(gene.coef,
                       freq,
                       noise.var,
                       Dskim = 0,
                       ...){
  gene.var <- calc.gene.var(gene.coef = gene.coef,
                            freq      = freq,
                            Dskim     = Dskim)
  herit <- gene.var/(gene.var + noise.var)
  return(herit)
}



##' Calculate the total genetic variance
##'
##' The genetic variance is calculated based on the genetic parameters
##' \eqn{\beta_{Gw_{t}}}{\beta_Gwt}, t in (1, 2, \eqn{\cdots}{...}, 8).
##'
##' as described in the publications by Kao and Cockerham:\cr \cr
##'
##' Kao CH, Zeng ZB. Modeling epistasis of quantitative trait loci using Cockerham's model. Genetics. 2002 Mar 1;160(3):1243-61. \cr
##' \url{http://www.genetics.org/content/160/3/1243.short} \cr \cr
##'
##' Cockerham CC, Weir BS. Quadratic analyses of reciprocal crosses. Biometrics. 1977 Mar 1:187-203. \cr
##' \url{http://www.jstor.org/stable/2529312}
##' @title Genetic Variance
##' @export
##' @param gene.coef a list with the coefficients of genetic effects.
##'                  Each element includes 8 genetic parameters (regression
##'                  coefficient)
##'                  \eqn{\beta_{Gw_{t}}}{\beta_Gwt}, t in (1, 2, \eqn{\cdots}{...}, 8)
##' @param freq a dataframe with the allele frequencies.
##' @param Dskim the coefficient of linkage disequilibrium. Default is 0 (no LD).
##' @param ... not used.
##' @return genetic variance.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de}
##' @examples
##' calc.gene.var(gene.coefficients, allele.freq)
calc.gene.var <- function(gene.coef,
                          freq,
                          Dskim = 0,
                          ...){
  gene.var <- 0
  for(i in 1: length(gene.coef)){
    epi <- gene.coef[[paste0("epi.par", i)]]
    fA <- freq[freq[,"SNP"] == epi$SNPA, "major.frequency"]
    fa <- freq[freq[,"SNP"] == epi$SNPA, "minor.frequency"]
    fB <- freq[freq[,"SNP"] == epi$SNPB, "major.frequency"]
    fb <- freq[freq[,"SNP"] == epi$SNPB, "minor.frequency"]

    taA <- fa - fA
    taB <- fb - fB
    piA <- fa * fA
    piB <- fb * fB
    a1 <- epi$additiveA
    d1 <- epi$dominanceA
    a2 <- epi$additiveB
    d2 <- epi$dominanceB
    iaa <- epi$additive_additive
    iad <- epi$additive_dominance
    ida <- epi$dominance_additive
    idd <- epi$dominance_dominance

    gene.var <- gene.var +
                2 * piA *a1^2 + piA * (1 + taA^2) * d1^2 + 2 * piB * a2^2 + piB * (1 + taB^2) * d2^2 +
                2 * (piA * taB^2 + piB * taA^2 + 2 * piA * piB - taA * taB * Dskim) * iaa^2 +
                0.25 * ((1-2*piA)-taB^2*(taA*taB + 4*Dskim)^2)*iad^2 +
                0.25 * ((1-2*piB)-taA^2*(taA*taB + 4*Dskim)^2)*ida^2 +
                0.0625 * (1-(taA*taB + 4*Dskim)^4)*idd^2 +
                4*piA*taA*a1*d1 + 4*Dskim*a1*a2 + 4*taB*Dskim*a1*d2 - 4*taB*piA*a1*iaa -
                2*(piA*taB^2 + 4*Dskim^2)*a1*iad - 2*(2*taB*piA*taA + (1-2*taA^2)*Dskim)*a1*ida -
                2*(taA*piA*taB^2 - (taA^2-4*piA)*taB*Dskim - 4*taA*Dskim^2)*a1*idd +
                4*taA*Dskim*d1*a2 - 4*Dskim*(taA*taB - 2*Dskim)*d1*d2 - 4*piA*(taA*taB+2*Dskim)*d1*iaa -
                2*taB*piA*(taA*taB + 4*Dskim)*d1*iad - 2*(taB*piA*(taA^2 + 1) + taA^3*Dskim)*d1*ida -
                (piA*taB^2*(1+taA^2) - 2*taA^3*taB*Dskim - 4*taA^2*Dskim^2)*d1*idd +
                4*taB*piB*a2*d2 - 4*taA*piB*a2*iaa -
                2*(2*taA*piB*taB + (1-2*taB^2)*Dskim)*a2*iad - 2*(piB*taA^2 + 4*Dskim^2)*a2*ida -
                2*(2*taB*piB*taA^2 - (taB^2 - 4*piB)*taA*Dskim - 4*taB*Dskim^2)*a2*idd -
                4*piB*(taA*taB + 2*Dskim)*d2*iaa -
                2*(taA*piB*(taB^2 + 1) + taB^3*Dskim)*d2*iad - 2*taA*piB*(taA*taB + 4*Dskim)*d2*ida -
                (piB*taA^2*(1 + taB^2) - 2*taB^3*taA*Dskim - 4*taB^2*Dskim^2)*d2*idd +
                2*(taB*(piA + 2*piB*taA^2) + taA*(1-3*taB^2)*Dskim - 4*taB*Dskim^2)*iaa*iad +
                2*(taA*(piB + 2*piA*taB^2) + taB*(1-3*taA^2)*Dskim - 4*taA*Dskim^2)*iaa*ida +
                0.5*(taA*taB + 2*Dskim)*(1 - (taA*taB + 4*Dskim)^2)*iaa*idd +
                0.5*(taA*taB + 2*Dskim - taA*taB*(taA*taB + 4*Dskim)^2)*iad*ida +
                0.25*(taA - taB*(taA*taB + 4*Dskim)^3)*iad*idd + 0.25*(taB - taA*(taA*taB + 4*Dskim)^3)*ida*idd
  }
  return(gene.var)
}


##' Give suggestion on the parameter setting of noise variance according to the expected heritability.
##' @title Suggestion noise
##' @export
##' @param gene.coef a list including the coefficients of genetic effects.
##'                  Each element includes 8 genetic parameters (regression
##'                  coefficient)
##'                  \eqn{\beta_{Gw_{t}}}, t in (1, 2, \eqn{\cdots}{...}, 8)
##' @param freq a dataframe with the allele frequencies.
##' @param Dskim the coefficient of linkage disequilibrium. Default is 0 (no LD).
##' @param heritability expected heritability.
##' @param ... not used.
##' @return variance of noise to generate the random noise.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de}
##' @examples
##' get.noise.var(gene.coefficients, allele.freq, 0.5)
get.noise.var <- function(gene.coef,
                          freq,
                          heritability,
                          Dskim = 0,
                          ...){
  gene.var <- calc.gene.var(gene.coef = gene.coef,
                            freq      = freq,
                            Dskim     = Dskim)
  noise.var <- gene.var * (1 - heritability) / heritability
  return(noise.var)
}



