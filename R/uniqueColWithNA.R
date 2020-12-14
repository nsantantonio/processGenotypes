#' uniqueColWithNA function
#'
#' function to (do something)
#'
#' @param X [value]
#' @param corThreshold [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
uniqueColWithNA <- function(X, corThreshold = 1){
	# Xmiss <- matrix(is.na(X), nrow = nrow(X), ncol = ncol(X))
	Xmiss <- is.na(X)
	nMiss <- colSums(Xmiss)

	suppressWarnings(corX <- cor(X, use = "pairwise.complete.obs"))
	corX <- corX >= corThreshold
	
	corXnoVar <- which(is.na(corX), arr.ind = TRUE)
	if (nrow(corXnoVar) > 0) {
		sameNoVar <- apply(corXnoVar, 1, function(x) all(X[, x[1]] == X[, x[2]], na.rm = TRUE))
		corX[corXnoVar] <- sameNoVar
	}
	
	corInd <- which(corX, arr.ind = TRUE)
	sameInd <- unique(t(apply(corInd[corInd[,1] != corInd[,2],], 1, sort)))
	if (nrow(sameInd) > 0) {
		mostMiss <- apply(sameInd, 1, function(x) x[which.max(nMiss[x])])
		return(!1:ncol(X) %in% unique(mostMiss))
	} else {
		return(rep(TRUE, ncol(X)))
	}	
}
