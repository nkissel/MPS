######################################################################################################
# naRemover, fills in NAs, likely only has use or applicability for full.select method
######################################################################################################
#' Remove NAs to MPS Graph Matrix
#'
#' Internal function that removes \code{NA}s to matrix that holds MPS graph data
#' @param myMat matrix of MPS data
#' @keywords matrix
#' @examples mat_input <- c(
#'   "Var1", NA, NA, NA, "Var2", NA,
#'   "Var3", NA, "Var4", "Var5", "Var4", "Var5"
#' )
#' mat1 <- matrix(mat_input, nrow = 4)
#' naRemover(mat1)
#' @export
naRemover <- function(myMat) {
  y <- dim(myMat)[2]
  for (i in 1:(y - 1)) {
    tv <- ifelse(is.na(myMat[, i]), "REP", myMat[, i])
    myRle <- rle(tv)
    toRep <- na.omit(ifelse(
      (myRle$values != "REP") &
        (c(myRle$values[-1], "HOLDER") != "REP"), myRle$lengths,
      ifelse(
        (myRle$values != "REP") &
          (c(myRle$values[-1], "HOLDER") == "REP"),
        c(myRle$lengths[-1], 1) + myRle$lengths, NA
      )
    ))
    myMat[, i] <- rep(myRle$values[myRle$values != "REP"], toRep)
  }
  return(myMat)
}

######################################################################################################
# naAdder, adds NAs, likely only has use or applicability for full.select method
######################################################################################################
#' Add NAs to MPS Graph Matrix
#'
#' Internal function that adds \code{NA}s to matrix that holds MPS graph data
#' @param myMat matrix of MPS data
#' @keywords matrix
#' @examples mat_input <- c(
#'   "Var1", "Var1", "Var1", "Var1", "Var2", "Var2",
#'   "Var3", "Var3", "Var4", "Var5", "Var4", "Var5"
#' )
#' mat1 <- matrix(mat_input, nrow = 4)
#' naAdder(mat1)
#' @export
naAdder <- function(myMat) {
  x <- dim(myMat)[1]
  y <- dim(myMat)[2]
  sticker <- rep(FALSE, x)
  if (y > 1) {
    for (i in 1:(y - 1)) {
      myRle <- rle(myMat[, i])
      toNA <- cumsum(myRle$lengths) + 1
      toNA <- c(1, toNA[-length(toNA)])
      indNA <- rep(FALSE, length(myMat[, i]))
      indNA[toNA] <- TRUE
      sticker <- sticker | (indNA)
      myMat[!sticker, i] <- NA
    }
  }
  return(myMat)
}
