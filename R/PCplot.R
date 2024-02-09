#' PCplot function
#'
#' function to (do something)
#'
#' @param X [value]
#' @param PC [value]. Default is NULL
#' @param xflip [value]. Default is FALSE
#' @param yflip [value]. Default is FALSE
#' @param pdfFileName [value]. Default is NULL
#' @param title [value]. Default is NULL
#' @param height [value]. Default is 7
#' @param width [value]. Default is 7
#' @param widthRatio [value]. Default is TRUE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
PCplot <- function(X, PC = NULL, xflip = FALSE, yflip = FALSE, pdfFileName = NULL, title = NULL, height = 7, width = 7, widthRatio = TRUE, ...){
	udv <- svd(X)
	if(is.null(PC)) PC <- 1:2
 	if(xflip) udv$u[, PC[1]] <- udv$u[, PC[1]] * -1
    if(yflip) udv$u[, PC[2]] <- udv$u[, PC[2]] * -1

    if(widthRatio) width <- round(udv$d[PC[1]] / udv$d[PC[2]] * height,2)
	# r <- cor(udv$u[, PC[1]] * udv$d[PC[1]], udv$u[, PC[2]] * udv$d[PC[2]])
	# print(r)
	txtplot(udv$u[, PC[1]] * udv$d[PC[1]], udv$u[, PC[2]] * udv$d[PC[2]])
	plot(udv$u[, PC[1]] * udv$d[PC[1]], udv$u[, PC[2]] * udv$d[PC[2]], 
    	    xlab = paste0("PC", PC[1], " ", round({udv$d^2}[1] / sum(udv$d^2) * 100, 2), "%"), 
        	ylab = paste0("PC", PC[2], " ",  round({udv$d^2}[2] / sum(udv$d^2) * 100, 2), "%"),
        	main = paste0("", title), ...)
    	# legend("topleft", legend = paste0("r = ", round(r, 2)))
}
