#' vanRaden function
#'
#' function to (do something)
#'
#' @param M [value]
#' @param coding [value]. Default is NULL
#' @param MARGIN [value]. Default is 2
#' @param center [value]. Default is TRUE
#' @param directCalcEpi [value]. Default is FALSE 
#' @param addDiag [value]. Default is 0.01
#' @param returnVar [value]. Default is FALSE
#' @param returnMAF [value]. Default is FALSE
#' @param returnCentered [value]. Default is FALSE
#' @param returnImputed [value]. Default is FALSE
#' @param checkPhase [value]. Default is FALSE
#' @param imp [value]. Default is "mean"
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
vanRaden <- function(M, coding = NULL, MARGIN = 2, center = TRUE, directCalcEpi = FALSE , addDiag = 0.01, returnVar = FALSE, returnMAF = FALSE, returnCentered = FALSE, returnImputed = FALSE, checkPhase = FALSE, imp = "mean"){
	mode <- function(x) {
		xmode <- table(x)
		as.numeric(names(xmode)[which.max(xmode)])
	}

	if(class(M) == "data.frame") M <- as.matrix(M)

	if(is.null(coding) || !coding %in% c("01", "101", "012")) {
		codes <- unique(c(M))
		if (all(c(-1, 1) %in% codes) | all(codes %in% c(-1, 0, 1))) {
			coding <- "101" 
		} else if (all(c(0, 2) %in% codes) | all(codes %in% c(0, 1, 2))) {
			coding <- "012"
		} else if(all(c(0, 1) %in% codes) | all(codes %in% c(0, 1))) {
			coding <- "01"
		} else {
			stop("Please provide a marker matrix of -101, 012, or 01 coding!")
		}
	} 

	if(coding %in% c("101", "012")) center <- TRUE

	if(MARGIN == 1) M <- t(M)

	n <- nrow(M)
	m <- ncol(M)

	getP <- function(M, code){
		if(code == "01"){
			p <- colMeans(M, na.rm = TRUE)
		} else {
			if(code == "101") p <- colMeans((M + 1) / 2, na.rm = TRUE) else p <- colMeans(M / 2, na.rm = TRUE) 
		}
		p
	}

	p <- getP(M, coding)

	if(checkPhase){
		wrongPhase <- p > 0.5
		if(any(wrongPhase)){
			cat("Warning! there are some alleles that have not been phased by MAF!\nAlleles will be phased by MAF...\n")
			if(coding == "01") {
				M[, wrongPhase] <- (M[, wrongPhase] - 1) * -1
			} else if(coding == "101"){
				M[, wrongPhase] <- (M[, wrongPhase]) * -1
			} else {
				M[, wrongPhase] <- (M[, wrongPhase] - 1) * -1 + 1
			}
			p <- getP(M, coding)
		}
	}
	
	alVar <- c(2 * crossprod(p, (1 - p)))

	NAs <- is.na(M)
	if(any(NAs)){
		cat("Marker matrix contains missing values! Imputing with population", imp, "\n")
		# if(imp == "mode" & !checkPhase) stop("you must specify 'checkPhase = TRUE' to use the population mode imputation")
		if(imp == "mean") {
			M <- apply(M, 2, function(x) {x[is.na(x)] <-  mean(x, na.rm = TRUE); x}) 
		} else {
			 M <- apply(M, 2, function(x) {x[is.na(x)] <- -1; x})
		}
	}

	
	if(center) Mcentered <- scale(M, scale = FALSE) else Mcentered <- M

	MMt <- tcrossprod(Mcentered)
	G <- MMt / alVar 
	if(coding == "01") G <- 2 * G

	if(addDiag){
		G <- G + diag(addDiag, n)
	}

	outList <- list(G = G)
	if(returnImputed) outList <- c(outList, list(M = M))
	if(returnCentered) outList <- c(outList, list(Z = Mcentered))
	if(returnVar) outList <- c(outList, alVars)
	if(returnMAF) outList <- c(outList, alFreqs)
	if(length(outList) == 1) return(G) else return(outList)
}
