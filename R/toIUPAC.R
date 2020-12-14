#' toIUPAC function
#'
#' function to (do something)
#'
#' @param X [value]
#' @param NAsymbol [value]. Default is "---"
#' @param checkX [value]. Default is TRUE
#' @param includePM [value]. Default is TRUE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
toIUPAC <- function(X, NAsymbol = "---", checkX = TRUE, includePM = TRUE){
	sortNucleo <- function(x) sapply(strsplit(x, ""), function(xx) paste0(sort(xx), collapse = ""))

	if(class(X) == "data.frame") {
		returnDf <- TRUE
		X <- as.matrix(X)
	} else {
		returnDf <- FALSE
	}

	iupac <- c(AA = "A", CC = "C", GG = "G", TT = "T", AC = "M", AG = "R", AT = "W", CG = "S", CT = "Y", GT = "K", CA = "M", GA = "R", TA = "W", GC = "S", TC = "Y", TG = "K")
	if(includePM){	
		mp <- expand.grid(c("+", "-"), c("A", "C", "G", "T", "+", "-"))
		ns <- apply(mp, 1, paste, collapse = "")
		pm <- c(ns, ns)
		names(pm) <- c(ns, apply(mp[, 2:1], 1, paste, collapse = ""))
		pm[pm == "-+"] <- "+-"
		pm[pm == "++"] <- "+"
		pm[pm == "--"] <- "-"
		iupac <- c(iupac, pm[!duplicated(names(pm))])
	}

	if(checkX){
		X[X == NAsymbol] <- NA
		Xcalls <- apply(X, 2, nchar)
		isTwo <- Xcalls != 2
		isTwo[is.na(isTwo)] <- FALSE
		if(any(isTwo)) stop("I cant deal with more than two alleles per individual, please fix me!")
	}
	
	for(i in names(iupac)){
		isi <- X == i
		cat("Percent of scores with value", i, ":", sum(isi, na.rm = TRUE) / length(c(X)), "\n")
		X[isi] <- iupac[i]	
	}
	
	if(returnDf) X <- data.frame(X)
	X
}
