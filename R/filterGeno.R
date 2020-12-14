#' filterGeno function
#'
#' function to (do something)
#'
#' @param M [value]
#' @param MARGIN [value]. Default is 2
#' @param MAF [value]. Default is 0.01
#' @param MISS [value]. Default is 0.5
#' @param HET [value]. Default is 0.05
#' @param returnMatrix [value]. Default is TRUE
#' @param returnStats [value]. Default is FALSE
#' @param returnImputed [value]. Default is FALSE
#' @param rmDupl [value]. Default is FALSE
#' @param rmDuplBy [value]. Default is NULL
#' @param maxGdiff [value]. Default is 0
#' @param maxiter [value]. Default is 10
#' @param checkAllMiss [value]. Default is FALSE
#' @param returnBool [value]. Default is FALSE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
filterGeno <- function(M, MARGIN = 2, MAF = 0.01, MISS = 0.5, HET = 0.05, returnMatrix = TRUE, returnStats = FALSE, returnImputed = FALSE, rmDupl = FALSE, rmDuplBy = NULL, maxGdiff = 0, maxiter = 10, checkAllMiss = FALSE, returnBool = FALSE){

	# checkinput <- function(X) {
	# 	allMinDuplBy <- all(colnames(M) %in% c(rmDuplBy))
	# 	allDuplByinM <- all(c(rmDuplBy) %in% colnames(M))
	# }
	checkAllNA <- function(X){
		isnaX <- is.na(X)
		allNA1 <- rowSums(isnaX) / ncol(X) == 1 
		allNA2 <- colSums(isnaX) / nrow(X) == 1 
		if(any(allNA2)) cat(sum(allNA2), " markers contain all missing data!\n")
		if(any(allNA1)) cat(sum(allNA1), " genotypes contain all missing data!\n")
		return(X[!allNA1, !allNA2])
	}
	
	dimFilter <- function(X, margin, het = HET, maf = MAF, miss = MISS){
		isnaX <- is.na(X)
		if(margin == 1) pMiss <- rowMeans(isnaX) else pMiss <- colMeans(isnaX)
		# pHet <- apply(X, margin, function(x) length(which(x == 0)) / sum(!is.na(x)))
		# pHet <- apply(X, margin, function(x) sum(x == 0) / sum(!is.na(x)))
		if(margin == 1) pHet <- rowSums(X == 0, na.rm = TRUE) / rowSums(!isnaX) else pHet <- colSums(X == 0, na.rm = TRUE) / colSums(!isnaX)

		percMissOk <- pMiss <= miss
		if(any(is.na(pMiss))) percMissOk[is.na(pMiss)] <- FALSE

		percHetOk <- pHet <= het
		if(any(is.na(pHet))) percHetOk[is.na(pHet)] <- FALSE
		
		if(margin == 1)	cat(sum(!percHetOk), " lines too heterozygous!    ", sum(!percMissOk), " lines missing too many genotype calls!\n" )

		if(maf < 1){
			# if(maf < 1) af <- 0.5 + apply(X, margin, function(x) if(any(!is.na(x))) sum(x, na.rm = TRUE) / (2 * sum(!is.na(x))) else 0) 
			# af <- 0.5 + apply(X, margin, function(x) if(any(!is.na(x))) sum(x, na.rm = TRUE) / (2 * sum(!is.na(x))) else 0) 
			# af <- 0.5 + colSums(X, na.rm = TRUE) / (2 * colSums(!isnaX))
			af <- 0.5 * (1 + colMeans(X, na.rm = TRUE))
			af[is.na(af)] <- 0

			wrongPhase <- af > 0.5
			if(any(wrongPhase)) af[wrongPhase] <- 1 - af[wrongPhase]	
			MAFOk <- af >= maf
			return(MAFOk & percMissOk & percHetOk )
		}

		return(percMissOk & percHetOk)
	}

	imputeMean <- function(x){
		mu <- mean(x, na.rm = TRUE)
		x[which(is.na(x))] <- mu
		return(x)
	}
	
	iterFilter <- function(filMGlist){
		filG <- dimFilter(M[, filMGlist[["filM"]]], margin = 1, het = gHET, miss = gMISS, maf = 1)
		filM <- dimFilter(M[filMGlist[["filG"]], ], margin = 2, het = mHET, miss = mMISS)
		return(list(filG = filG, filM = filM))
	}
	classM <- class(M)
	if(!classM %in% c("data.frame", "matrix")) stop("M must be a data.frame or matrix!")
	if(classM %in% "data.frame") M <- as.matrix(M)

	if (MARGIN == 1) M <- t(M)
	if(checkAllMiss) M <- checkAllNA(M)

	if (length(HET) > 1){
		gHET <- HET[2]
		mHET <- HET[1]
	} else {
		mHET <- gHET <- HET
	}

	if (length(MISS) > 1){
		gMISS <- MISS[2]
		mMISS <- MISS[1]
	} else {
		mMISS <- gMISS <- MISS
	}

	n <- nrow(M)
	p <- ncol(M)
	filM <- dimFilter(M, margin = 2, het = mHET, miss = mMISS)
	filG <- dimFilter(M, margin = 1, het = gHET, miss = gMISS, maf = 1)
	
	if(any(filG)){
		filMGlisti <- NULL
		filMGlist <- list(filG = filG, filM = filM)
		cat("iteration = 0", "\n") 
		print(unlist(lapply(filMGlist,sum)))
		i <- 1
		continue <- TRUE
		while(continue & i <= maxiter){
			cat("iteration = ", i, "\n")
			i <- i + 1
			filMGlisti <- filMGlist 
			filMGlist <- iterFilter(filMGlist)
			print(unlist(lapply(filMGlist, sum)))
			Gdiff <- sum(filMGlisti[["filG"]]) - sum(filMGlist[["filG"]])
			continue <- (Gdiff < 0 | Gdiff > maxGdiff) 
		}
		cat("final dimensions: \n")
		print(unlist(lapply(filMGlist, sum)))
		filM <- filMGlist[["filM"]]
		filG <- filMGlist[["filG"]]
	}

	if(returnBool & !returnMatrix) return(filMGlist)
	M <- M[filG, filM]

	if (rmDupl){
		if(is.null(rmDuplBy)) {
			keepUnique <- uniqueColWithNA(M)
		} else {
			uniqueList <- list()
			for(i in 1:length(rmDuplBy)) {
				msubset <- colnames(M) %in% rmDuplBy[[i]]
				nsub <- sum(msubset)
				cat("split:\t", names(rmDuplBy)[i], "no. markers:\t", nsub, "\t")
				uniqueList[[i]] <- uniqueColWithNA(M[, msubset])
				psub <- sum(uniqueList[[i]])
				cat("no. markers removed:\t", nsub - psub, "\t", ", no. markers remaining:\t", psub, "\n")
			}
			keepUnique <- unlist(uniqueList)
		}
		M <- M[, keepUnique]
		cat("Total after removing duplicates: ", sum(keepUnique), "\n")
		
		filMGlist[["filM"]] <- 1:p %in% (1:p)[filMGlist[["filM"]]][keepUnique]
	}
	

	percMiss <- colSums(is.na(M)) / nrow(M)
	percHet <- apply(M, 2, function(x) length(which(x == 0))  / sum(!is.na(x)))

	alleleFreq <- 0.5 + apply(M, 2, function(x) sum(x, na.rm = TRUE) / (2 * sum(!is.na(x))))
	
	wrongPhase <- alleleFreq > 0.5
	wrongPhase[is.na(wrongPhase)] <- FALSE
	if(any(wrongPhase)) alleleFreq[wrongPhase] <- 1 - alleleFreq[wrongPhase]
	M[,wrongPhase] <- M[,wrongPhase] * -1

	if (returnImputed & any(is.na(M))) M <- apply(M, 2, imputeMean)

	if (MARGIN == 1) M <- t(M)
	if(classM %in% "data.frame") M <- as.data.frame(M)

	outList <- list(M = M)
	if(returnBool) outList <- c(outList, filMGlist)
	if (returnStats) outList <- c(outList, list(maf = alleleFreq, percMiss = percMiss, percHet = percHet))
	if (length(outList) == 1) return(outList[[1]]) else return(outList)
}
