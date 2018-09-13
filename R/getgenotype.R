##' Read genotype data.
##'
##' If it is plink file format (.bed, .bim, .fam),
##' make sure that plink has already been installed in the system
##' @title Read genotype data based on the file format
##' @export
##' @param fname a string specifying the file to read genotype information from
##' @param verbose when set show the commands that are to be called through \code{system}.
##' @param run when set (default) execute the \code{system} calls.
##' @param cleanup when set (default) remove intermediate files before returning.
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
##' @param ... not used.
##' @return a dataframe of genotype data: columns are the SNPs; rows are indviduals.
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de} and
##'         Benno Pütz \email{puetz@@psych.mpg.de}
##' @examples
##' ## "snp.head" genotype file: rows are individuals and columns are SNPs
##' # get full path of example file
##' fgeno.path <- system.file("extdata", "10SNP.txt", package="SimPhe")
##'
##' geno <- read.geno(fgeno.path, ftype = "snp.head")
##' head(geno)
##'
##' ## "plink" genotype file: 1).map and .ped; 2).bed, .fam and .bim
##' # get directory of plink example file
##' fpath <- strsplit(fgeno.path, "10SNP.txt")
##'
##' #### Note: before run this example, specify your installation path of plink ####
##' # geno <- read.geno(paste0(fpath, "bdemo"), ftype = "plink", plink.path = "user's plink path")
read.geno <- function(fname      = NULL,
                      verbose    = getOption('verbose'),
                      run        = TRUE,
                      cleanup    = TRUE,
                      ftype      = c("ind.head", "plink", "snp.head"),
                      plink.path = NULL,
                      ...){
  if (is.null(fname)){
      stop('please specify the name of the genotype file to read')
  }

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

  geno <- switch(ftype,
                 ind.head = {
                   return(t(read.table(file             = fname,
                                       header           = TRUE,
                                       row.names        = 1,
                                       stringsAsFactors = FALSE)))
                 },
                 plink = {
                   ## prepare
                   outfile1 <- paste0(fname,"_recode")
                   # detect current operating system
                   current.sys <- Sys.info()[['sysname']]

                   ## get plink path
                   if(!is.character(plink.path)){
                     if(current.sys %in% "Windows"){
                       plink.path <- system("where plink", intern = TRUE)
                     }else{
                       plink.path <- system("which plink", intern = TRUE)
                     }
                   }



                   # check plink.path
                   if(is.character(plink.path)){
                     if(file.exists(paste0(fname, ".bed"))){
                       cmd1 <- paste(plink.path, "--bfile", fname, "--noweb --recode A-transpose --out", outfile1)
                     }else
                     {
                       cmd1 <- paste(plink.path, "--file", fname, "--noweb --recode A-transpose --out", outfile1)
                     }
                   }else{
                     stop('please specify the installion path of plink')
                   }

                   if(verbose) cat(cmd1, file = stderr())
                   if(run) {
                     system(cmd1)
                     if(cleanup)
                       unlink(paste0(outfile1, '.log'))
                     outfile1 <- paste0(outfile1, '.traw') # add suffix used by plink (output file name of the --recode A-transpose)
                   }

                   if(current.sys %in% "Windows"){
                     ftraw <- read.table(file = outfile1,
                                         header = TRUE,
                                         row.names = 2,
                                         stringsAsFactors = FALSE)
                     if(cleanup)
                       unlink(outfile1)
                     return(t(ftraw[,c(-1:-5)]))
                   }else{
                     ##
                     outfile2 <- paste0(fname,".geno")
                     cmd2 <- paste("awk '{", '$1="";$3="";$4="";$5="";$6="";print $0',"}'",outfile1,
                                   "|",
                                   "sed -e 's/    / /g' -e 's/^ //' >", outfile2)
                     if(verbose) cat(cmd2, file = stderr())
                     if(run) {
                       system(cmd2)
                       ##
                       return(t(read.table(file = outfile2,
                                           header = TRUE,
                                           row.names = 1,
                                           stringsAsFactors = FALSE)))
                       if(cleanup){
                         unlink(outfile1)
                         unlink(outfile2)
                       }
                     }
                   }},
                 snp.head = {
                   return(read.table(file             = fname,
                                     header           = TRUE,
                                     row.names        = 1,
                                     stringsAsFactors = FALSE))
                 }
  )
  return(geno)
}
