##' Main process of simulation for phenotypes
##'
##' further discussion on pattern
##' @title Simulation for phenotypes (SimPhe main process)
##' @export
##' @param sim.pars a prepared list containing the parameters settings for simulation or a file of parameters settings.
##'                    Please set your own parameters
##'                    following the same structure as the object \code{genepars} or as the file
##'                    \code{simupars.txt}
##'                    (you could find example file \cr
##'                    \code{system.file("extdata","simupars.txt",package="SimPhe")}).\cr\cr
##'                    To specify heritability, there
##'                    are two ways: one is to set
##'                    heritability in parameter file or in the prepared list object which will futher pass to \code{sim.pars}.
##'                    Another way is to set \code{noise.var}
##'                    by using function \code{\link{get.noise.var}}
##'                    given specify \code{heritability}. \cr\cr
##'
##'                    List format: please follow the example of the object \code{genepars}.
##'                    The meaning of each elment in the list is similar like the file format description below.\cr\cr
##'
##'                    File format: please follow the example of the \code{simupars.txt} file
##'                    found in the inst/extdata/ directory of the package (run \cr
##'                    \code{system.file("extdata", "simupars.txt", package="SimPhe")} \cr
##'                    to to get the path to the file),
##'                    blank lines are ignored.
##'                    The file consists of three or four blocks for each phenotype (the number of blocks depends on user): mean, main and epistasis, sometimes heritability. Each block
##'                    is started by a line of the form '[blockname]' followed by the parameter setting for the block, e.g. for first phenotype,
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
##'                     \item{SNPA}{first SNP}
##'                     \item{SNPB}{second SNP}
##'                     \item{additive_additive}{coefficient for additive-additive interaction}
##'                     \item{additive_dominance}{coefficient for additive-dominance interaction}
##'                     \item{dominance_additive}{coefficient for dominance-additive interaction}
##'                     \item{dominance_dominance}{coefficient for dominance-dominance interaction}
##'                     }}
##'   \item{[P1heritability]}{\describe{
##'                         \item{heritability}{expected heritability}
##'                         }}
##' }
##'                    Similar meanings for "[P2mean]", "[P2main]", "[P2epistasis]", and so on.
##'
##' @param fgeno file to read genotype information from or pre-read data.frame
##'                    with that information (matching the output format of
##'                    \code{\link{read.geno}}).
##' @param ftype genotype file format, it accepts three options:
##' \describe{
##' \item{"plink":}{plink format (\code{.bed}, \code{.bim}, \code{.fam} or \code{.map}, \code{.ped});}
##' \item{"ind.head":}{columns are the indviduals and lines are SNPs;}
##' \item{"snp.head":}{columns are SNPs and lines are indviduals.}}
##'                    For "plink", \code{fgeno} needs to be given
##'                    without suffix and \code{plink.path} may need to be
##'                    assigned by the user because \code{plink} will be run
##'                    from within \pkg{SimPhe}. More detail see \code{plink.path}.
##'                    For the other options, \code{fgeno} should be the full
##'                    name (with suffix and path if necessary) of the
##'                    genotype file.
##'                    Of course, this does not apply if \code{fgeno} is
##'                    provided as a data frame.
##' @param plink.path  path of plink executable. Only needed when the \code{ftype}
##' is "\code{plink}". Default is NULL. The function will detect the plink path with
##' \code{system("where plink")} for Windows users and \code{system("which plink")} for
##' Linux and MacOS users. But there is no garantee that the commands work on all devices.
##' If the path cannot be determined or the executable cannot be called from \code{read.geno},
##' then users have to try other formats.
##' @param fwrite logical. Write out file (simulated data) or not. If TRUE (default),
##'                    simulated phenotypes will be written, respectively.
##' @param fphename filename of the phnotype(s). Default is "simu.pheno".
##' @param fusepar filename of the setting for simulation (for recording).
##'                    Default is "usedpars.txt".
##' @param Dskim the coefficient of linkage disequilibrium.
##'                    Default is 0 (no LD).
##' @param seed an integer used for set.seed(). Default is NA.
##' @param noise.var variance for random noise. Default is 1. Note that this
##'                   is overridden by the heritability setting in the
##'                   simulation parameter file. If heritability is
##'                   given in parameter file then \code{noise.var} will not work.
##' @param pattern ignore pattern for detecting the phenotype index from the parameter names.
##'                Default is "[[:alpha:]]+" which means letters.
##' @param genetic.model a string show the genetic model to use for simulation.
##'                   Default is "epistasis".
##' @param ... not used.
##' @return a data.frame with the simulated phenotype(s) where the column(s)
##'                   refer to different phenotype(s) and rows to individuals.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de} and
##'                   Benno Pütz \email{puetz@@psych.mpg.de}
##' @examples
##' #### file path of example
##' # simulation parameters:
##' fpar.path <- system.file("extdata", "simupars.txt", package="SimPhe")
##'
##' # genotype file: rows are individuals and columns are SNPs
##' fgeno.path <- system.file("extdata", "10SNP.txt", package="SimPhe")
##'
##'
##' #### instead of a parameter file, prepared list like genepars also works
##' genepars
##'
##'
##' #### simulate phenotype(s)
##' phe <- sim.phe(sim.pars = fpar.path, fgeno = fgeno.path, ftype = "snp.head", fwrite = FALSE)
##' # or
##' phe <- sim.phe(sim.pars = genepars, fgeno = fgeno.path, ftype = "snp.head", fwrite = FALSE)
##'
##' # the simulated phenotype(s)
##' str(phe)
##' head(phe)
sim.phe <- function(sim.pars = NULL,
                    fgeno         = NULL,
                    ftype         = c("ind.head", "plink", "snp.head"),
                    fwrite        = TRUE,
                    fphename      = "simu.pheno",
                    fusepar       = "usedpars.txt",
                    seed          = NA,
                    Dskim         = 0,
                    noise.var     = 1,
                    pattern       = "[[:alpha:]]+",
                    plink.path    = system("which plink"),
                    genetic.model = "epistasis",
                    ...){
  if(!is.na(seed))
    set.seed(seed)

  ftype <- tryCatch(match.arg(ftype),
                    error = function(e) {
                      print(e)
                      print(class(e))
                      frmls <- eval(formals(sys.function(1))[['ftype']])
                      warning("cannot handle ftype '", ftype,
                              "', possible options are:\n\t",
                              paste0(frmls, collapse  = ', '),
                              "\n\nWill use default")
                      return(frmls[1])
                    }
  )
  if(is.list(sim.pars)){
    genetic.pars <- sim.pars
  }else{
    genetic.pars <- read.simu.pars(file = sim.pars)    # test for NULL there
  }


  nphe <- length(unique(na.omit(as.numeric(unlist(strsplit(names(genetic.pars),
                                                           split = pattern))))))
  if (!is.data.frame(fgeno)) {
    geno <- read.geno(fname      = fgeno,               # test for NULL there
                      ftype      = ftype,
                      plink.path = plink.path)

  } else {
    geno <- fgeno
  }
  check.snp.par(genetic.pars = genetic.pars,
                nphe         = nphe)


  phe <- list()
  for (i in 1:nphe) {
    # calculate frequency
    genetic.pars[[paste0("P", i, "frequency")]] <- get.freq(geno     = geno,
                                                            epi.pars = specify.pars(genetic.pars = genetic.pars,
                                                                                    effect.type  = "epistasis",
                                                                                    phe.index    = i))
    gene.coef <- get.gene.coef(main.pars = specify.pars(genetic.pars = genetic.pars,
                                                        effect.type  = "main",
                                                        phe.index    = i),
                               epi.pars  = specify.pars(genetic.pars = genetic.pars,
                                                        effect.type  = "epistasis",
                                                        phe.index    = i),
                               model     = genetic.model)

    if(!is.null(genetic.pars[[paste0("P", i, "heritability")]])){
      exp.noise.var <- get.noise.var(gene.coef    = gene.coef,
                                     freq         = genetic.pars[[paste0("P", i, "frequency")]],
                                     Dskim        = Dskim,
                                     heritability = as.numeric(genetic.pars[[paste0("P", i, "heritability")]]))
    } else {
      exp.noise.var <- noise.var
      # calculate heritability
      genetic.pars[[paste0("P", i, "heritability")]] <- calc.herit(gene.coef = gene.coef,
                                                                   freq      = genetic.pars[[paste0("P", i, "frequency")]],
                                                                   noise.var = exp.noise.var)
    }

    phe[[paste0("p", i)]] <- as.numeric(genetic.pars[[paste0("P", i, "mean")]]) +
      gene.effect(geno      = geno,
                  gene.coef = gene.coef,
                  model     = genetic.model) +
      rnorm(nrow(geno),
            sd = sqrt(exp.noise.var)) ## now automatically but need to think more
    rownames(phe[[paste0("p", i)]]) <- rownames(geno)
  }
  pars.writer(genetic.pars = genetic.pars,
              fname        = fusepar)
  phe <- list2frame(phe)
  if(fwrite)
    phe.writer(phe   = phe,
               fname = fphename)
  return(phe)
}






##' Get genetic effect for each individual based on the genotype.
##' @title Get genetic effect
##' @export
##' @param geno a data.frame or a matrix containing the genotype information.
##' @param gene.coef a list with the coefficients of genetic effects.
##' @param model a string specifying the genetic model to use for the simulation.
##'               Default is "epistasis".
##' @param ... not used.
##' @return a data.frame including genetic effects.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de} and
##'         Benno Pütz \email{puetz@@psych.mpg.de}
##' @examples
##' # genotype file: rows are individuals and columns are SNPs
##' fgeno.path <- system.file("extdata", "10SNP.txt", package="SimPhe")
##'
##' # get genotype
##' geno <- read.geno(fgeno.path, ftype = "snp.head")
##'
##' # take a look at geno and gene.coef
##' geno
##' gene.coefficients
##'
##' # get gene effects
##' gene.effect(geno, gene.coefficients)
gene.effect <- function(geno,
                        gene.coef,
                        model = c("epistasis"),
                        ...){
  gene.effect <- data.frame(rep(0, nrow(geno)))
  model <- tryCatch(match.arg(model),
                    error = function(e) {
                      print(e)
                      print(class(e))
                      frmls <- eval(formals(sys.function(1))[['model']])
                      warning("cannot handle model '", model,
                              "', possible options are:\n\t",
                              paste0(frmls, collapse  = ', '),
                              "\n\nWill use default")
                      return(frmls[1])
                    }
  )
  switch(model,
         epistasis = {
           for (i in 1:nrow(geno)){ # each individual
             for (j in 1:length(gene.coef)){ # each interaction
               eff.scale <- as.matrix(genetic.scale(SNPA = geno[i, gene.coef[[paste0("epi.par", j)]]$SNPA],
                                                    SNPB = geno[i, gene.coef[[paste0("epi.par", j)]]$SNPB]))
               gene.effect[i,] <- (gene.effect[i,] +
                                     (as.matrix(gene.coef[[paste0("epi.par", j)]][c(-1, -2)]) %*% eff.scale))
             }
           }
         },
         ## handle any unknown model specifications here
         {
           stop('should have never come here ...')
         }
  )
  return(gene.effect)
}






##' Get the coefficients of genetic effectsre.
##' @title Get the coefficients of genetic effect
##' @export
##' @param main.pars a data.frame or a matrix containing the parameters for the  main effect:
##'                  additive and dominace.
##' @param epi.pars a data.frame or a matrix containing the parameters for the epistatic effect:
##' additive  \eqn{\times}{*} additive,
##' additive  \eqn{\times}{*} dominance,
##' dominance \eqn{\times}{*} additive,
##' dominance \eqn{\times}{*} dominance.
##' @param model a string show the genetic model to use for simulation. Default is "epistasis"
##' @param ... not used.
##' @return a list with the coefficients of genetic effects.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de} and
##'                   Benno Pütz \email{puetz@@psych.mpg.de}
##' @examples
##' # take a look at the settings of coefficients for main effects
##' maineff.pars
##'
##' # take a look at the settings of coefficients for interactive effects
##' epistasis.pars
##'
##' # get a vector of gene coefficients
##' get.gene.coef(maineff.pars, epistasis.pars)
get.gene.coef <- function(main.pars,
                          epi.pars,
                          model = c("epistasis"),
                          ...){
  coef <- list()

  model <- tryCatch(match.arg(model),
                    error = function(e) {
                      print(e)
                      print(class(e))
                      frmls <- eval(formals(sys.function(1))[['model']])
                      warning("cannot handle model '", model,
                              "', possible options are:\n\t",
                              paste0(frmls, collapse  = ', '),
                              "\n\nWill use default")
                      return(frmls[1])
                    }
  )

  switch(model,
         epistasis = {
           for(j in 1:nrow(epi.pars)){
             coef[[paste0("epi.par", j)]] <- data.frame(SNPA             = epi.pars[j, "SNPA"],
                                                        SNPB             = epi.pars[j, "SNPB"],
                                                        stringsAsFactors = FALSE)

             gene.eff <- unlist(c(main.pars[main.pars$SNP == epi.pars[j, "SNPA"], -1],
                                  main.pars[main.pars$SNP == epi.pars[j, "SNPB"], -1],
                                  epi.pars[j, 3:6]))
             coef[[paste0("epi.par", j)]] <- cbind(coef[[paste0("epi.par", j)]],
                                                   matrix(gene.eff, nrow = 1))
             names(gene.eff)[1:2] <- paste0(names(gene.eff)[1:2], "A")
             names(gene.eff)[3:4] <- paste0(names(gene.eff)[3:4], "B")
             names(coef[[paste0("epi.par", j)]]) <- c("SNPA", "SNPB", names(gene.eff))
           }
         },
         ## unknown model specifications have been handles above
         {
           stop('should have never come here ...')
         }
  )
  return(coef)
}


##' Get the parameters of main/epistatic effects per phenotype.
##' @title Get the parameters of main/epistatic effects per phenotype
##' @export
##' @param genetic.pars a data.frame or a matrix containing the parameters
##'                     information for main effect: additive and dominance.
##' @param effect.type a string naming the type of the genetic effects
##'                     (accepts either "main" or "epistasis").
##' @param phe.index a integer indicating the phenotype. Default is 1.
##' @param ... not used.
##' @return a data.frame or a matrix containing the parameters information for epistatic effect:
##' additive  \eqn{\times}{*} additive,
##' additive  \eqn{\times}{*} dominance,
##' dominance \eqn{\times}{*} additive,
##' dominance \eqn{\times}{*} dominance.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de}
##' @examples
##' # get parameters of coefficients for main effects
##' specify.pars(genepars, effect.type = "main")
##'
##' # get parameters of coefficients for interactive effects
##' specify.pars(genepars, effect.type = "epistasis")
specify.pars <- function(genetic.pars,
                         effect.type = c("main", "epistasis"),
                         phe.index   = 1,
                         ...){
  effect.type <- tryCatch(match.arg(effect.type),
                          error = function(e) {
                            print(e)
                            print(class(e))
                            frmls <- eval(formals(sys.function(1))[['effect.type']])
                            warning("cannot handle effect.type '", effect.type,
                                    "', possible options are:\n\t",
                                    paste0(frmls, collapse  = ', '),
                                    "\n\nWill use default")
                            return(frmls[1])
                          }
  )

  effect.type.options <- eval(formals()[['effect.type']])
  if(effect.type %in% effect.type.options)
    return(genetic.pars[[paste0("P", phe.index, effect.type)]])
  else
    stop("please specify the effect type, should be one of\n\t",
         paste(effect.type.options, collapse = ', '))
}

