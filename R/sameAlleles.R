#' sameAlleles function
#'
#' function to (do something)
#'
#' @param x [value]
#' @param y [value]
#' @param returnTables [value]. Default is FALSE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
sameAlleles <- function(x, y, returnTables = FALSE){
	xall <- lapply(strsplit(x, "\\/"), sort)
	yall <- lapply(strsplit(y, "\\/"), sort)
	
	xlen <- sapply(xall, length)
	ylen <- sapply(yall, length)

	xcall <- sapply(xall, paste, collapse = "/")
	ycall <- sapply(yall, paste, collapse = "/")
	allSame <- xcall == ycall
	cat(sum(allSame) / length(allSame), "% same allele calls\n")

	sum(xlen != ylen)
	
	alleleConflict <- list()
	alleleConflictCounter <- 0
	xallNotInyall <- list()
	xallNotInyallCounter <- 0
	yallNotInxall <- list()
	yallNotInxallCounter <- 0

	notSame <- which(!allSame)
	
	for(i in notSame){
		if(xlen[i] == ylen[i]){
			alleleConflictCounter <- alleleConflictCounter + 1
			alleleConflict[[alleleConflictCounter]] <- c(index = i, allelesInCommon = sum(xall[[i]] %in% yall[[i]]))
		} else if (xlen[i] < ylen[i]) {
			xallNotInyallCounter <- xallNotInyallCounter + 1
			xiny <- xall[[i]] %in% yall[[i]]
			xallNotInyall[[xallNotInyallCounter]] <- c(index = i, allelesInCommon = sum(xiny), allxIny = all(xiny))
		} else {
			yallNotInxallCounter <- yallNotInxallCounter + 1
			yinx <- yall[[i]] %in% xall[[i]]
			yallNotInxall[[yallNotInxallCounter]] <- c(index = i, allelesInCommon = sum(yinx), allyInx = all(yinx))
		}
	}

	hardConflict <- data.frame(do.call(rbind, alleleConflict))
	xinyDf <- data.frame(do.call(rbind, xallNotInyall))
	yinxDf <- data.frame(do.call(rbind, yallNotInxall))

	yNotInx <- yinxDf[yinxDf[["allyInx"]] == 0, ]
	xNotIny <- xinyDf[xinyDf[["allxIny"]] == 0, ]

	thirdAlleleInx <- yinxDf[yinxDf[["allelesInCommon"]] == 2, ]
	thirdAlleleIny <- xinyDf[xinyDf[["allelesInCommon"]] == 2, ]

	list(allSame = which(allSame), yNotInx = yNotInx$index, xNotIny = xNotIny$index, thirdAlleleInx = thirdAlleleInx$index, thirdAlleleIny = thirdAlleleIny$index)
}
