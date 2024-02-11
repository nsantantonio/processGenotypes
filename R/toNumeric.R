#' toNumeric function
#'
#' function to (do something)
#'
#' @param x [value]
#' @param check [value]. Default is FALSE
#' @param hetAltAlleleToHomoMinor [value]. Default is TRUE
#' @param returnAf [value]. Default is TRUE
#' @param includePM [value]. Default is TRUE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
toNumeric <- function(x, check = FALSE, hetAltAlleleToHomoMinor = TRUE, returnAf = TRUE, includePM = TRUE, inclHasHet = TRUE){
	iupac <- c(AA = "A", CC = "C", GG = "G", TT = "T", AC = "M", AG = "R", AT = "W", CG = "S", CT = "Y", GT = "K") 
	homo <- c("A", "C", "G", "T")
	hetero <- c("K", "Y", "W", "S", "R", "M")
	if(inclHasHet) {
		hetero <- c("H", hetero)
		iupac <- c(iupac, H = "H")
	}
	if(includePM) {
		mpHet <- c(apply(expand.grid(c("+", "-"), c("A", "C", "G", "T")), 1, paste, collapse = ""), "+-")
		pm <- mpHet
		names(pm) <- pm
		pm["+"] <- "+"
		pm["-"] <- "-"
		iupac <- c(iupac, pm)
		homo <- c(homo, "+", "-")
		hetero <- c(hetero, mpHet)
	}
	
	xtab <- table(x)
	xhomo <- xtab[names(xtab) %in% homo]
	xhetero <- xtab[names(xtab) %in% hetero]

	alleleCounter <- xhomo * 2
	for(i in names(xhetero)){
		for(j in strsplit(names(iupac)[iupac %in% i], "")[[1]]) alleleCounter[j] <- alleleCounter[j] + xhetero[i]
	}
	alleleCounter <- sort(alleleCounter, decreasing = TRUE)
	af <- alleleCounter / sum(alleleCounter)
	major <- names(alleleCounter)[1]
	minor <- names(alleleCounter)[2]

	if(length(xhomo) > 2 | length(xhetero) > 1){
		# cat("More than two alleles! Setting third (smallest) allele to minor...\n ")
		altAlleles <- sort(names(xhomo)[names(xhomo) != major])
		if (hetAltAlleleToHomoMinor) x[x %in% iupac[paste(altAlleles, collapse = "")]] <- minor #set heterozygous non-major to homozygous minor
		x[x %in% altAlleles] <- minor # set homozygous alternate to minor	
	}
	x[x == major] <- -1
	x[x == minor] <- 1
	x[x %in% names(xhetero)] <- 0
	if (check) {
		if (length(xhomo) > 1 & length(xhetero) == 1) {
			if(iupac[paste(sort(names(xhomo)), collapse = "")] != xhetero) stop("Wrong heterozygous call! Need to investigate further...")
		} else if(length(xhomo) > 1) stop("Too many heterozygous type calls! I cant handle 3 or more nucleotide snps. Please fix me!")
	}

	if(returnAf) return(list(x = x, alleleFreq = af, genoTable = xtab)) else return(x)  
}
