##' Write out the information about the parameters for simulation.
##' @title Write current parameters to file
##' @export
##' @param genetic.pars list of simulation parameters. One element per block of
##'                     \code{file}. Each element is a dataframe with SNP names
##'                     and model parameters.
##' @param fname filename of the setting for simulation (for recording).
##' @param ... not used.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de}
##' @examples
##' pars.writer(genepars)
pars.writer <- function(genetic.pars,
                        fname = "usedpars.txt",
                        ...){
  if(file.exists(fname))
    file.remove(fname) # if parmeters file exist, delete it and rewrite it
  heads <- paste("[", names(genetic.pars), "]", sep = "")
  for (i in 1:length(genetic.pars)) {
    write.table(heads[i],
                file      = fname,
                append    = TRUE,
                quote     = FALSE,
                row.names = FALSE,
                col.names = FALSE)
    write.table(genetic.pars[[i]],
                file      = fname,
                append    = TRUE,
                quote     = FALSE,
                row.names = FALSE,
                col.names = FALSE)
  }
}



##' Write out the simulated phenotypes.
##' @title Write phenotypes
##' @export
##' @param phe a data.frame or a matrix of simulated phenotypes.
##' Each column is a phenotype.
##' @param onefile whether to create just one file for all phenotypes (default) or
##'                one per phenotype
##' @param fname filename of the phnotype(s).
##' @param ... not used.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de} and
##'         Benno Pütz \email{puetz@@psych.mpg.de}
##' @examples
##' phe <- matrix(rnorm(1000), ncol = 2)
##' colnames(phe) <- c("p1", "p2")
##' phe.writer(phe)
phe.writer <- function(phe,
                       onefile = TRUE,
                       fname   = "simu.pheno",
                       ...)
{
  if (has.rn <- !is.null(rownames(phe))){
    phe <- cbind(ID = rownames(phe),
                 phe)
  }
  if(onefile){
    write.table(phe,
                file      = fname,
                row.names = FALSE,
                col.names = TRUE,
                quote     = FALSE)
  }else{
    cnp <- colnames(phe)
    for (i in 1:(ncol(phe) - has.rn)){
      col.set <- which(cnp %in% c('ID', cnp[i+has.rn]))
      write.table(phe[, col.set, drop = FALSE],
                  file      = fname,
                  row.names = FALSE,
                  col.names = TRUE,
                  quote     = FALSE)
    }
  }
}


##' Convert list to data.frame.
##' @title Convert list to data.frame
##' @export
##' @param x a list. In this package, it is used as a list includes simulated
##' phenotypes. One element per block.
##' Each element is a dataframe with SNP names and individuals.
##' @param ... not used.
##' @return a data.frame.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de}
##' @examples
##' x <- list(test1=matrix(rnorm(1000), ncol=2), test2=matrix(rnorm(1000), ncol=2))
##' str(x)
##' x.new <- list2frame(x)
##' str(x)
list2frame <- function(x,
                       ...){
  namecol <- names(x)
  x <- data.frame(x)
  colnames(x) <- namecol
  return(x)
}
