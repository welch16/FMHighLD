##' @useDynLib FMHighLD
##' @importFrom Rcpp evalCpp
##' @exportPattern "^[[:alpha:]]+"
##' @import GenomicRanges
NULL


#' @rdname FMDataSet
#' @export
setClass("FMDataSet",
         contains = "GRanges",
         representation = representation(
           annot_conf = "data.frame",
           response = "data.frame",
           LD_mat = "matrix"))

setValidity("FMDataSet", function(object) {
  # if (! ("counts" %in% assayNames(object)) )
  #   return( "the assays slot must contain a matrix named 'counts'" )
  # if ( !is.numeric( counts(object) ) )
  #   return( "the count data is not numeric" )
  # if ( any( is.na( counts(object) ) ) )
  #   return( "NA values are not allowed in the count matrix" )
  # if ( !is.integer( counts(object) ) )
  #   return( "the count data is not in integer mode" )
  # if ( any( counts(object) < 0 ) )
  #   return( "the count data contains negative values" )
  # 
  # design <- design(object)
  # # 'design' is either a formula or matrix
  # stopifnot(is(design, "formula") | is(design, "matrix"))
  # 
  # if (is(design, "formula")) {
  #   designVars <- all.vars(design)
  #   if (!all(designVars %in% names(colData(object)))) {
  #     return("all variables in design formula must be columns in colData")
  #   }
  #   designVarsClass <- sapply(designVars, function(v) class(colData(object)[[v]]))
  #   if (any(designVarsClass == "character")) {
  #     return("variables in design formula are character vectors.
  #            convert these columns of colData(object) to factors before including in the design formula")
  #   }
  #   designFactors <- designVars[designVarsClass == "factor"]
  #   # levels would duplicate after make.names()
  #   if (any(sapply(designFactors,function(v) {
  #     factor.lvls <- levels(colData(object)[[v]])
  #     factor.nms <- make.names(factor.lvls)
  #     any(duplicated(factor.nms))
  #   }))) {
  #     return("levels of factors in the design have non-unique level names after make.names() is applied.
  #            best to only uobject letters and numbers for levels of factors in the design")
  #   }
  #   # levels contain characters other than letters, numbers, and underscore
  #   if (any(sapply(designFactors,function(v) {
  #     factor.lvls <- levels(colData(object)[[v]])
  #     any(!grepl("^[A-Za-z0-9_.]+$",factor.lvls))
  #   }))) {
  #     # just a warning for now
  #     message("  Note: levels of factors in the design contain characters other than
  #             letters, numbers, '_' and '.'. It is recommended (but not required) to use
  #             only letters, numbers, and delimiters '_' or '.', as these are safe characters
  #             for column names in R. [This is a message, not an warning or error]")
  #   }
  #   } else if (is(design, "matrix")) {
  #     # TODO add some more tests for if 'design' is matrix
  #     stopifnot(nrow(design) == ncol(object))
  #   }
  # 
  TRUE
    })

##' FMDataSet object and constructors
##'
##' \code{FMDataSet} is a a subclass of \code{GenomicRanges}, used to store the 
##' input values, intermediate calculations and results of a multi-response fine-mapping analysis 
##' of SNPs in High LD.
##' 
##' @param snp_data a data.frame with at least the columns SNP, seqnames and position. 
##' The remaining columns are going to be considered as annotations.
##' @param annot_conf a data.frame with the configuration of the annotation data as fixed or 
##' random effects.
##' @param response a data.frame with three columns SNP, pheno and the value of an association 
##' statistic between both. For example, for a continuous phenotype we would use the wald statistic.
##' @param LD_matrix a matrix with the correlation between all the SNPs in snp_data.
##' 
##' @return A FMDataSet object.
##' 
##' @aliases FMDataSet FMDataSet-class 
##' @docType class
##' @rdname FMDataSet
##' @importFrom utils packageVersion
##' @export
##' @examples 
#' snps =c("a","b","c","d","e","f")
#' nsnps = length(snps)
#' nassoc = 20
#' snp_data = data.frame(SNP = snps,
#'                   seqnames = rep("chrA",nsnps),
#'                   position = sample(1e6,nsnps),
#'                   annot1 = factor(sample(c(TRUE,FALSE),nsnps,replace = TRUE )),
#'                   annot2 = rnorm(nsnps))
#' annot_conf = data.frame(annot = c("annot1","annot2"),
#'                         fixed = c(TRUE,TRUE),random = c(TRUE,FALSE))
#' response = data.frame(SNP = sample(snps,nassoc,replace = TRUE),
#'                      pheno = sample(c("I","II","III"),nassoc,replace = TRUE),
#'                      z = rnorm(nassoc,sd = 30))
#' LD_matrix = matrix(runif(nsnps^2),nrow = nsnps)
#' LD_matrix = .5 * (LD_matrix + t(LD_matrix))
#'
#' fm = FMDataSet(snp_data,annot_conf,response,LD_matrix)
#' fm
FMDataSet <- function(snp_data,annot_conf,response,LD_matrix)
{
  stopifnot(c("seqnames","position","SNP") %in% colnames(snp_data))

  ranges = with(snp_data,
                GRanges(seqnames = seqnames,
                        ranges = IRanges(
                          start = position,
                          width = 1
                        )))
  snps = snp_data$SNP
  names(ranges) = snps
  annot_data = dplyr::select(snp_data,-seqnames,-position,-SNP)
  mcols(ranges) = as(annot_data,"DataFrame")
  
  stopifnot("annot" %in% colnames(annot_conf))
  coldata = column_to_rownames(annot_conf,"annot")

  new("FMDataSet",ranges,annot_conf = coldata,response = response,LD_mat = LD_matrix)
}
