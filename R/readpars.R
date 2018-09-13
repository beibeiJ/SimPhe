##' Extract (sub)strings matching regex pattern.
##'
##' Extract the substrings of \code{x} that match the regex pattern \code{pattern}.
##' The pattern may contain groups (enclosed in parentheses) which will result
##' in further substrings extracted.
##'
##' Derived from the help on \code{\link[base]{regmatches}}, take a look at the help there
##' @title Extract (sub)strings matching regex pattern
##' @export
##' @param x a character vector.
##' @param pattern regular expression to be found in \code{x}.
##' @param ... not used.
##' @return a character (string) matrix where the first column contains the global
##' match for the pattern, each pair of '()' will result in another column with
##' the respective match.
##' @author Benno Pütz \email{puetz@@psych.mpg.de}
##' @examples
##' s <- "Test: A1 BC23 DEF456"
##' pattern = "([[:alpha:]]+)([[:digit:]]+)"
##' regextract(s, pattern)
##'
##' # equivalent to this example from the help page for grep()
##' lapply(regmatches(s, gregexpr(pattern, s)), function(e) regmatches(e, regexec(pattern, e)))
regextract <- function(x,
                       pattern,
                       ...){
  n <- length(x)
  return(matrix(unlist(sapply(regmatches(x,
                                         gregexpr(pattern, x)),
                              function(e) regmatches(e, regexec(pattern, e)))),
                nrow  = n,
                byrow = TRUE))
}




##' Read file specifying the simulation parameters.
##'
##' File format: please follow the example of the \code{simupars.txt} file
##' found in the inst/extdata/ directory of the package (run \cr
##' \code{system.file("extdata", "simupars.txt", package="SimPhe")} \cr
##' to see the full path),
##' blank lines are ignored.
##' The file consists of three or four blocks for each phenotype:
##' mean, main, epistasis, and (optionally) heritability.
##' Each block is started by a line of the form '[blockname]' followed by the
##' parameters for the block, e.g., for thefirst phenotype,
##' \describe{
##'   \item{[P1mean]}{\describe{
##'                 \item{mean}{\eqn{\beta_{0}}{\beta_0}: coefficient parameter of "basic" genetic effects in
##'                            \eqn{G_{ij}=\beta_{0}+\sum_{t=1}^{8}\beta_{G_{w_{t}}}w_{tij}}{Gij = \beta_0 + \sum(\beta_Gwt*wtij), t in (1, 2, ..., 8)}. }
##'                 }}
##'   \item{[P1main]}{\describe{
##'                 \item{SNP}{SNP name}
##'                 \item{additive}{coefficient of additive effect}
##'                 \item{dominance}{coefficient of dominance effect}
##'                 }}
##'   \item{[P1epistasis]}{\describe{
##'                     \item{SNPA}{first SNP name}
##'                     \item{SNPB}{second SNP name}
##'                     \item{additive_additive}{coefficient for additive-additive interaction}
##'                     \item{additive_dominance}{coefficient for additive-dominance interaction}
##'                     \item{dominance_additive}{coefficient for dominance-additive interaction}
##'                     \item{dominance_dominance}{coefficient for dominance-dominance interaction}
##'                     }}
##'   \item{[P1heritability]}{\describe{
##'                         \item{heritability}{expected heritability}
##'                         }}
##' }
##' For each block the expected columns and the respective menaings are given.
##' Similar blocks need to be provided for the other phenotype(s): "[P2mean]",
##' "[P2main]", "[P2epistasis]", and so on.
##' @title Read parameters
##' @export
##' @param file the file of parameters settings
##' @param ... not used
##' @return a list of simulation parameters. One element per block of \code{file}.
##' Each element is a dataframe with SNP names and model parameters
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de} and
##'         Benno Pütz \email{puetz@@psych.mpg.de}
##' @examples
##' # simulation parameters:
##' fpar.path <- system.file("extdata", "simupars.txt", package="SimPhe")
##'
##'
##' # pars <- read.simu.pars(fpar.path)
read.simu.pars <- function(file = NULL,
                           ...){
  if (is.null(file)){
    ## once we have a functional package we can use system.file()
    ## to automatically read the accompanying simupars.txt file -
    ## for now there is no secure way to determine its path
    stop("please specify the parameter file to read")
  }
  #detect operating system
  current.sys <- Sys.info()[['sysname']]
  if(current.sys %in% "Windows"){
    file.new <- pipe(paste('sed -e "/^$/d"', file))
  }else{
    ## initial parsing to get block structure information
    file.new <- pipe(paste('sed -e /^$/d', file))
  }
  all.lines <- readLines(file.new)
  close(file.new)
  # all.lines <- readLines(pipe(paste('sed -e /^$/d', file)))

  block.starts <- grep('^\\[', all.lines)          # all block header lines

  block.nrows <- diff(c(block.starts,
                        length(all.lines)+1)) - 2  # one for header and one for
                                                   # [block] header

  ## second parsing of the individual blocks
  pars <- list()                          # start with empty list
  for (i in 1:length(block.starts)) {
    block.name <- regextract(all.lines[block.starts[i]], '[[:alnum:]]+')
    pars[[block.name]] <- read.table(textConnection(all.lines),
                                     header           = TRUE,
                                     stringsAsFactors = FALSE,
                                     blank.lines.skip = TRUE,
                                     skip             = block.starts[i],
                                     nrows            = block.nrows[i])
  }
  return(pars)
}


##' Check the number of the SNPs set in the simulation parameters.
##'
##' The number of SNPs involved in main effects should be the same as the number
##' of SNPs involved in epistasis.
##' @title Check the number of the SNPs involved in epistasis and main effects
##' @export
##' @param genetic.pars a data.frame or a matrix containing the parameter
##'                     information for the main effect: additive and dominance.
##' @param nphe number of phenotypes.
##' @param ... not used.
##' @return NULL --- will stop if test fails.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de}
##' @examples
##' check.snp.par(genepars, nphe = 2)
check.snp.par <- function(genetic.pars,
                          nphe,
                          ...){
  for (i in 1:nphe) {
    nsnp.main <- length(genetic.pars[[paste0("P", i, "main")]][, "SNP"])
    nsnp.epi <- length(genetic.pars[[paste0("P", i, "epistasis")]][, "SNPA"]) +
                length(genetic.pars[[paste0("P", i, "epistasis")]][, "SNPB"])
    if(nsnp.epi != nsnp.main)
      stop("Please check the parameter file and make sure the number of SNPs involved in main effects
           is same as the number of SNPs involved in epistasis!")
  }
}

# gp <- function(i, name){
#   genetic.pars[[paste0("P", i, name)]]
# }
