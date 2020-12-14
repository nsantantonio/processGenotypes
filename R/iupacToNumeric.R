#' iupacToNumeric function
#'
#' function to (do something)
#'
#' @param X [value]
#' @param NAsymbol [value]. Default is "N"
#' @param snpMargin [value]. Default is 2
#' @param checkX [value]. Default is TRUE
#' @param returnMAF [value]. Default is FALSE
#' @param includePM [value]. Default is TRUE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
iupacToNumeric <- function(X, NAsymbol = "N", snpMargin = 2, checkX = TRUE, returnMAF = FALSE, includePM = TRUE){
	homo <- c("A", "C", "G", "T")
	hetero <- c("K", "Y", "W", "S", "R", "M")
	if(includePM) {
		mpHet <- c(apply(expand.grid(c("+", "-"), c("A", "C", "G", "T")), 1, paste, collapse = ""), "+-")
		homo <- c(homo, "+", "-")
		hetero <- c(hetero, mpHet)
	}

	if (class(X) == "data.frame") {
		returnDf <- TRUE
		X <- as.matrix(X)
	} else {
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
	if(returnDf) M <- data.frame(M) 
	if(returnMAF) return(list(X = M, alleleFreq = alleleFreq, genoTable = genoTable)) else return(M)
}
