##' Frequencies of SNPs.
##'
##' A dataset containing sample allele frequencies of SNPs.
##' allele.freq.
##'
##' @format A data frame with 6 rows and 3 columns(variables):
##' \describe{
##'   \item{SNP}{name of SNP}
##'   \item{major.frequency}{frequency of major allele}
##'   \item{minor.frequency}{frequency of minor allele}
##' }
##' @docType data
##' @keywords datasets
##' @name allele.freq
##' @usage data(allele.freq)
NULL



##' Coefficients of genetic effects.
##'
##' A dataset containing the regression coefficients of genetic effects.
##' gene.coefficients.
##'
##' @format A list with 3 elements:
##' \describe{
##' \item{epi.par1}{a data frame with 1 row and 10 variables:
##'                 \describe{
##'                   \item{SNPA}{first SNP}
##'                   \item{SNPB}{second SNP}
##'                   \item{additiveA}{coefficient for additive effect at locus A}
##'                   \item{dominanceA}{coefficient for dominance effect at locus A}
##'                   \item{additiveB}{coefficient for additive effect at locus B}
##'                   \item{dominanceB}{coefficient for dominance effect at locus B}
##'                   \item{additive_additive}{coefficient for additive-additive interaction}
##'                   \item{additive_dominance}{coefficient for additive-dominance interaction}
##'                   \item{dominance_additive}{coefficient for dominance-additive interaction}
##'                   \item{dominance_dominance}{coefficient for dominance-dominance interaction}
##'                   }}
##' \item{epi.par2}{a data frame with 1 row and 10 variables:
##'                 \describe{
##'                   \item{SNPA}{first SNP}
##'                   \item{SNPB}{second SNP}
##'                   \item{additiveA}{coefficient for additive effect at locus A}
##'                   \item{dominanceA}{coefficient for dominance effect at locus A}
##'                   \item{additiveB}{coefficient for additive effect at locus B}
##'                   \item{dominanceB}{coefficient for dominance effect at locus B}
##'                   \item{additive_additive}{coefficient for additive-additive interaction}
##'                   \item{additive_dominance}{coefficient for additive-dominance interaction}
##'                   \item{dominance_additive}{coefficient for dominance-additive interaction}
##'                   \item{dominance_dominance}{coefficient for dominance-dominance interaction}
##'                   }}
##' \item{epi.par3}{a data frame with 1 row and 10 variables:
##'                 \describe{
##'                   \item{SNPA}{first SNP}
##'                   \item{SNPB}{second SNP}
##'                   \item{additiveA}{coefficient for additive effect at locus A}
##'                   \item{dominanceA}{coefficient for dominance effect at locus A}
##'                   \item{additiveB}{coefficient for additive effect at locus B}
##'                   \item{dominanceB}{coefficient for dominance effect at locus B}
##'                   \item{additive_additive}{coefficient for additive-additive interaction}
##'                   \item{additive_dominance}{coefficient for additive-dominance interaction}
##'                   \item{dominance_additive}{coefficient for dominance-additive interaction}
##'                   \item{dominance_dominance}{coefficient for dominance-dominance interaction}
##'                   }}
##' }
##' @docType data
##' @keywords datasets
##' @name gene.coefficients
##' @usage data(gene.coefficients)
NULL




##' Parameter settings for simulation.
##'
##' A dataset containing the parameters for a simulation.
##' genepars.
##'
##' @format A list with 7 elements:
##' \describe{
##' \item{P1mean}{A data frame with 1 row and 1 variable:
##'               \describe{
##'                 \item{mean}{\eqn{\beta_{0}}: coefficient parameter of "basic" genetic effects in \eqn{G_{ij}=\beta_{0}+\sum_{t=1}^{8}\beta_{G_{w_{t}}}w_{tij}}. }
##'                 }}
##' \item{P1main}{A data frame with 6 rows and 3 variables:
##'               \describe{
##'                 \item{SNP}{SNP name}
##'                 \item{additive}{coefficient of additive effect}
##'                 \item{dominance}{coefficient of dominance effect}
##'                 }}
##' \item{P1epistasis}{A data frame with 3 rows and 6 variables:
##'                   \describe{
##'                     \item{SNPA}{first SNP}
##'                     \item{SNPB}{second SNP}
##'                     \item{additive_additive}{coefficient for additive-additive interaction}
##'                     \item{additive_dominance}{coefficient for additive-dominance interaction}
##'                     \item{dominance_additive}{coefficient for dominance-additive interaction}
##'                     \item{dominance_dominance}{coefficient for dominance-dominance interaction}
##'                     }}
##' \item{P1heritability}{A data frame with 1 row and 1 variable:
##'                      \describe{
##'                         \item{heritability}{expected heritability}
##'                         }}
##' \item{P2mean}{A data frame with 1 row and 1 variable:
##'               \describe{
##'                 \item{mean}{mean of genetic effect}
##'                 }}
##' \item{P2main}{A data frame with 6 rows and 3 variables:
##'               \describe{
##'                 \item{SNP}{SNP name}
##'                 \item{additive}{coefficient of additive effect}
##'                 \item{dominance}{coefficient of dominance effect}
##'                 }}
##' \item{P2epistasis}{A data frame with 3 rows and 6 variables:
##'                   \describe{
##'                     \item{SNPA}{first SNP}
##'                     \item{SNPB}{second SNP}
##'                     \item{additive_additive}{coefficient for additive-additive interaction}
##'                     \item{additive_dominance}{coefficient for additive-dominance interaction}
##'                     \item{dominance_additive}{coefficient for dominance-additive interaction}
##'                     \item{dominance_dominance}{coefficient for dominance-dominance interaction}
##'                     }}
##' }
##' @docType data
##' @keywords datasets
##' @name genepars
##' @usage data(genepars)
NULL


##' Parameter settings of main effects (additive effect and dominance).
##'
##' A dataset containing the parameter settings for main effects.
##' maineff.pars.
##'
##' @format A data frame with 6 rows and 3 columns(variables):
##' \describe{
##'   \item{SNP}{name of SNP}
##'   \item{additive}{coefficient of additive effect}
##'   \item{dominance}{coefficient of dominance effect}
##' }
##' @docType data
##' @keywords datasets
##' @name maineff.pars
##' @usage data(maineff.pars)
NULL


##' Parameter settings of epistatic effects.
##'
##' A dataset containing the parameter settings for epistatic effects.
##' epistasis.pars.
##'
##' @format A data frame with 6 rows and 3 columns(variables):
##' \describe{
##'   \item{SNPA}{first SNP}
##'   \item{SNPB}{second SNP}
##'   \item{additive_additive}{coefficient for additive-additive interaction}
##'   \item{additive_dominance}{coefficient for additive-dominance interaction}
##'   \item{dominance_additive}{coefficient for dominance-additive interaction}
##'   \item{dominance_dominance}{coefficient for dominance-dominance interaction}
##' }
##' @docType data
##' @keywords datasets
##' @name epistasis.pars
##' @usage data(epistasis.pars)
NULL



