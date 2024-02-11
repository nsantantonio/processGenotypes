#' iupacToNumeric function
#'
#' function to (do something)
#'
#' @param X data.frame or matrix of IUPAC genotype calls 
#' @param NAsymbol character for missing call symbol. Default is "N"
#' @param snpMargin integer. margin for snps. the other margin should be for individuals. 1 for rows, 2 for columns. Default is 2
#' @param checkX logical. Should X be checked for > 2 alleles per site? Default is TRUE
#' @param returnMAF logical. Should tables of allele and genotype frequencies be returned? Default is FALSE
#' @param includePM logical. Should '+' and '-' symbols be processed? Default is TRUE. 
#' @return data.frame or matrix (depending on input) of genotype scores as numeric values.
#' @details This funciton expects single character IUPAC scores for each genotype call. FOr example, an 'AA' call would be 'A', and an 'AT' call would be 'W'
#' @examples none
#' @export
iupacToNumeric <- function(X, NAsymbol = "N", snpMargin = 2, checkX = TRUE, returnMAF = FALSE, includePM = TRUE, inclHasHet = TRUE){
	homo <- c("A", "C", "G", "T")
	hetero <- c("K", "Y", "W", "S", "R", "M")
	if(inclHasHet) {
		hetero <- c("H", hetero)
	}
	if(includePM) {
		mpHet <- c(apply(expand.grid(c("+", "-"), c("A", "C", "G", "T")), 1, paste, collapse = ""), "+-")
		homo <- c(homo, "+", "-")
		hetero <- c(hetero, mpHet)
	}

	if ("data.frame" %in% class(X)) {
		returnDf <- TRUE
		X <- as.matrix(X)
	} else {
		if(!"matrix" %in% class(X)) stop(paste0("expecting a data.frame or matrix of snp scores. Got class ", class(X)[[1]]))
		returnDf <- FALSE
	}

	X[X == NAsymbol] <- NA
	if (checkX) {
		vals <- apply(X, snpMargin, table)
		if(snpMargin == 1) names(vals) <- rownames(X)
		uniqVals <- unique(unlist(sapply(vals, names)))
		if(!all(uniqVals %in% homo | uniqVals %in% hetero)) stop("Something is wrong... there are symbols other than IUPAC nucleotides. I cant handle 3 way heterozygote snps or other weird symbols. Please fix the input or me!")
		nGeno <- sapply(vals, length)
		nGenoTab <- table(nGeno)
		if(any(nGeno > 3)){
			cat("These markers have more than 2 alleles. Setting third (smallest) allele to minor. Use 'returnMAF = TRUE' to return genotypic tables and allele frequencies. \n\n")
			print(names(vals[nGeno > 3]))
		}
	}

	M <- apply(X, snpMargin, toNumeric, returnAf = returnMAF)
	if(returnMAF) {
		alleleFreq <- lapply(M, "[[", "alleleFreq")
		genoTable <- lapply(M, "[[", "genoTable")		
		if(snpMargin == 1) M <- do.call(rbind, lapply(M, "[[", "x")) else M <- do.call(cbind, lapply(M, "[[", "x"))
		dim(M)
	} else {
		if(snpMargin == 1) M <- t(M)
	}

	class(M) <- "numeric"
	if(returnDf) M <- as.data.frame(M) 
	if(returnMAF) return(list(X = M, alleleFreq = alleleFreq, genoTable = genoTable)) else return(M)
}
