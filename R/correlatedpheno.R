##' Build correlated phenotypes
##' @title Build correlated phenotypes
##' @export
##' @param pheno a matrix or dataframe with the phenotypic information.
##' @param corMtr a correlation matrix.
##' @param sdMtr a matrix with the standard deviation,
##'     e.g., if the number of dimensions is 2, then it is
##'     \eqn{\left(\begin{array}{cc} \sigma_{1} & 0 \\ 0 & \sigma_{2} \end{array}\right)}.
##'     If it is NULL (default), generate it based on the data from \code{pheno}.
##' @param margin a vector giving the subscript which the function will be
##' applied over.
##' E.g., for a matrix 1 indicates rows, 2 (default) indicates columns.
##' Where \code{pheno} has named dimnames, it can be a character vector selecting
##' dimension names.
##' @param ... not used.
##' @return a matrix with correlated phenotypes.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de}
##' @examples
##' x1 <- rnorm(4000, mean = 5, sd = 10)
##' x2 <- rnorm(4000, mean = 10, sd = 30)
##' x <- matrix(cbind(x1, x2), ncol = 2)
##'
##' # test original correlation
##' cor.test(x[, 1], x[, 2])
##'
##' # correlation matrix
##' corM <- matrix(c(1, 0.6, 0.6, 1), ncol = 2)
##'
##' # standard deviation matrix
##' sdM <- matrix(c(10, 0, 0, 30), ncol = 2)
##'
##' # build correlation
##' x.new <- build.cor.phe(x, corM, sdM)
##'
##' # check mean and standard deviation of new data set
##' apply(x.new, 2, mean)
##' apply(x.new, 2, sd)
##'
##' # test correlation
##' cor.test(x.new[, 1], x.new[, 2])
build.cor.phe <- function(pheno,
                          corMtr,
                          sdMtr        = NULL,
                          margin       = 2,
                          ...){
  # mean.type <- match.arg(mean.type)
  if(!(margin %in% 1:2)){
    stop("invalid value for parameter 'margin'")
  }

  if(margin == 1){
    pheno <- t(pheno)
    margin <- 2
  }
  nphe <- ncol(pheno)
  if(is.null(sdMtr))
    sdMtr <- build.sd.matrix(x = pheno)

  u <- apply(pheno, margin, mean)

  L <- chol(corMtr)
  cor.phe.stan <- t(L) %*% t(scale(pheno))
  cor.phe.var <- sdMtr %*% cor.phe.stan
  cor.phe <- matrix(ncol = nphe, nrow = nrow(pheno))
  for (i in 1:nphe)
    cor.phe[,i] <- cor.phe.var[i,] + u[i]
  colnames(cor.phe) <- paste0("p", 1:nphe)
  return(cor.phe)
}


##' Build a matrix with standard deviation.
##' @title build a matrix with standard deviation
##' @export
##' @param x a matrix or dataframe.
##' @param margin an integer specifying the dimension that standard
##'               deviation will be computed for.
##'               1 indicates rows, 2 indicates columns. Default is 2.
##' @param ... not used
##' @return a matrix with the standard deviation.
##'     If the number of dimensions is 2, then it is
##'     \eqn{\left(\begin{array}{cc} \sigma_{1} & 0 \\ 0 & \sigma_{2} \end{array}\right)}.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de} and
##'         Benno Pütz \email{puetz@@psych.mpg.de}
##' @examples
##' x1 <- rnorm(4000, mean =  5, sd = 10)
##' x2 <- rnorm(4000, mean = 10, sd = 30)
##' x <- matrix(cbind(x1, x2), ncol = 2)
##' build.sd.matrix(x)
build.sd.matrix <- function(x,
                            margin = 2,
                            ...){
  # sdm <- matrix(rep(0,
  #                   nrow = ncol(x),
  #                   ncol = ncol(x))
  # diag(sdm) <- apply(x, margin, sd)
  sdm <- diag(apply(x, margin, sd))

  return(sdm)
}



