# imports stats graphics

#' @importFrom stats as.formula na.omit
#' @importFrom graphics lines par rect segments strheight strwidth text

####### THESE FUNCTIONS NEED TO BE REWRITTEN AS THEYRE VERY DIFFICULT TO PARSE

################################################################################
# Read Tree Label, likely only has use for build.tree function
################################################################################
#' @title Tree Processing Function
#'
#' @description This is an internal function used for building MPS graph.
#' @param var.mat.na matrix of MPS data
#' @keywords matrix
#' @export
read_tree_label <- function(var.mat.na) {
  # The inner workings of this method are actually somewhat complex bc use of RLE and the
  #   existence of so many NA. Rest assured I'll eventually get around to documenting
  #   all of this better, but for the time being, I won't add anything

  ##################################### data input and read ####################
  if (is.null(dim(var.mat.na))) {
    endPoint <- length(var.mat.na)
  } else {
    endPoint <- dim(var.mat.na)[2]
  }
  if (is.null(dim(var.mat.na))) {
    txt <- (var.mat.na[endPoint])
  } else {
    txt <- (var.mat.na[, endPoint])
  }
  if (is.null(dim(var.mat.na))) {
    v1 <- ifelse(is.na(var.mat.na[endPoint - 1]), FALSE, TRUE)
  } else {
    v1 <- ifelse(is.na(var.mat.na[, endPoint - 1]), FALSE, TRUE)
  }
  if (dim(var.mat.na)[2] > 1) {
    c <- 0
    same <- 0
    for (i in 1:length(v1)) {
      if (v1[i]) {
        c <- c + 1
      }
      same[i] <- c
    }
    rle1 <- rle(same)
    txt2 <- rep(NA, length(rle(same)$values))
    c <- 1
    for (i in 1:length(txt2)) {
      txt2[i] <- paste("(", paste(txt[c:(c + rle(same)$lengths[i] - 1)], collapse = ","), ")", sep = "")
      c <- c + rle(same)$lengths[i]
    }
    if (endPoint > 2) {
      for (h in (endPoint - 2):1) {
        if (is.null(dim(var.mat.na))) {
        	v2 <- ifelse(is.na(var.mat.na[h]), FALSE, TRUE)
        } else {
        	v2 <- ifelse(is.na(var.mat.na[, h]), FALSE, TRUE)
        }
        c <- 0
        same2 <- 0
        for (i in seq_along(v2)) {
          if (v2[i]) {
            c <- c + 1
          }
          same2[i] <- c
        }
        rle2 <- rle(same2)

        txt3 <- rep(NA, length(rle2$values))
        rle.summer <- (rle1$lengths[1])
        if (length(rle1$lengths) > 1) {
          for (i in 2:length(rle1$lengths)) {
            rle.summer[i] <- rle.summer[i - 1] + (rle1$lengths[i])
          }
        }

        adjuster <- 0
        c <- 0
        for (i in seq_along(txt3)) {
          txt3[i] <- paste("(", paste(txt2[(c + 1):((which((rle.summer - adjuster) == rle2$lengths[i])))], collapse = ","), ")", sep = "")
          c <- which((rle.summer - adjuster) == rle2$lengths[i])
          adjuster <- adjuster + rle2$lengths[i]
        }
        txt2 <- txt3
        rle1 <- rle2
      }
    } else {
      txt3 <- txt2
    }
  } else {
    txt3 <- paste(as.vector(var.mat.na))
  }
  final <- paste("(", paste(txt3, collapse = ","), ");", sep = "") # no-var at top
  return(final)
}


################################################################################
################################################################################
########################  Build Tree Graphics Function  ########################
################################################################################
################################################################################
#' @title Tree Graphic Function
#'
#' @description Builds MPS graphic in graphics window
#' @param var.mat.na matrix of MPS data
#' @param cex text size
#' @param alt integer that vertically staggers terminal nodes (max of 3)
#' @param shape1 either \code{'trad'} or \code{'rad'} for traditional trees and radial trees
#' @param shape2 either \code{'right'} or \code{'tri'} for rectangular and triangular;
#'  only matters for \code{shape1='trad'}
#' @param col node color
#' @keywords tree
#' @examples mat_input <- c(
#'   "Var1", "Var1", "Var1", "Var1", "Var2", "Var2",
#'   "Var3", "Var3", "Var4", "Var5", "Var4", "Var5"
#' )
#' mat1 <- matrix(mat_input, nrow = 4)
#' build.tree(mat1) # base example
#' build.tree(mat1, alt = 2) # different alt
#' build.tree(mat1, shape1 = "rad") # radial tree
#' build.tree(mat1, shape2 = "tri") # triangular tree (only for shape1)
#' build.tree(mat1, col = "red") # node color
#' @export
build.tree <- function(var.mat.na, cex = .9, alt = 1, shape1 = "trad", shape2 = "right", lwd = 1, col) {
  ############################ input info and stops ############################
  if (!is.numeric(alt)) alt <- 1
  if (shape1 == "rad") shape2 <- NULL
  if (shape1 != "trad" && shape1 != "rad") stop("invalid shape")
  if (shape2 != "right" && !is.null(shape2)) shape2 <- "tri" # only matters for shape1="right"
  if (missing(col)) col <- "lightgoldenrod"

  ############################ begin read and setup ############################
  var.mat.na <- naAdder(var.mat.na)
  final <- read_tree_label(var.mat.na)
  par(mar = c(0, 0, 0, 0))
  # plot(NA,NA,xlim=c(-1.1,1.1),ylim=c(-1.1,1.1), axes = TRUE)
  final <- gsub("_", "FILLER", final)
  final <- gsub("\\.", "FILLERPERIOD", final)
  str <- unlist(strsplit(final, "(?<=\\pP)(?=\\PP)|(?<=\\PP)(?=\\pP)", perl = T))
  maxDepth <- dim(var.mat.na)[2]
  x.depth <- rep(NA, length(na.omit(as.vector(var.mat.na))))
  current <- 0
  str <- str[-length(str)]
  shift.right <- 1
  shift <- 1

  ####################### Node-information matrix building #####################
  for (i in 1:length(str)) {
    shift.up.down <- 0
    temp.str <- unlist(strsplit(str[i], "(?<=\\pP)(?=\\pP)", perl = T))
    current.var <- sum(!is.na(x.depth)) + 1
    if (i == 1) {
      shift.up.down <- 1:length(unlist(strsplit(str[i], "(?<=\\pP)(?=\\pP)", perl = T)))
      x.depth[1:length(shift.up.down)] <- shift.up.down
      shift[1:length(shift.up.down)] <- 1
      current <- maxDepth + 1
    } else {
      if (("," %in% temp.str) & (length(temp.str) > 1)) {
        shift.right <- shift.right + 1
        shift.up.down <- length(which(!temp.str == ",")) / 2
        current <- current - shift.up.down
        x.depth[current.var:(current.var + (shift.up.down) - 1)] <- current:(maxDepth)
        shift[current.var:(current.var + (shift.up.down) - 1)] <- shift.right
        current <- (maxDepth)
      } else {
        if (!"," %in% temp.str) {
          x.depth[current.var] <- (maxDepth + 1)
          shift[current.var] <- shift.right
          current <- (maxDepth + 1)
        }
      }
    }
  }
  x.depth.true <- na.omit(x.depth)
  tot.mat <- (matrix(ncol = 6, nrow = length(x.depth.true)))
  tot.mat[, 1] <- c(1:length(x.depth.true))
  tot.mat[, 2] <- as.numeric(x.depth.true)
  tot.mat[, 3] <- as.numeric(shift)
  tot.mat[, 6] # pointer to father variable

  ######################### Terminal node coordinates ##########################
  if (shape1 == "rad") {
    r <- seq(from = 0, to = 1, len = maxDepth + 2 + alt - 1)
    r <- r[-1]
    theta <- seq(from = 0, to = 2 * pi, len = length(which(x.depth.true == maxDepth + 1)) + 1)
    theta <- theta[-length(theta)]
  }
  if (shape1 == "trad") {
    xTerm <- seq(-.95, .95, len = length(which(x.depth.true == maxDepth + 1)))
    if (dim(var.mat.na)[1] == 2 && dim(var.mat.na)[2] == 1) {
    	xTerm <- seq(-.4, .4, len = length(which(x.depth.true == maxDepth + 1)))
    }
    if (length(which(x.depth.true == maxDepth + 1)) == 1) xTerm <- 0
    y.heights <- seq(to = -1, from = 1, len = maxDepth + 2)
    y.heights <- y.heights[-length(y.heights)]
    if (alt == 2) {
      y.heights <- c(
        y.heights,
        y.heights[length(y.heights)] - strheight("1", "user") * 2 * cex
      )
    }
    if (alt == 3) {
      y.heights <- c(
        y.heights,
        y.heights[length(y.heights)] - strheight("1", "user") * 2 * cex,
        y.heights[length(y.heights)] - strheight("1", "user") * 4 * cex
      )
    }
  }

  ####################### Internal node coordinates ############################
  for (i in unique(rev(x.depth.true))) {
    if (i == maxDepth + 1) {
      if (shape1 == "trad") {
        tot.mat[which(tot.mat[, 2] == maxDepth + 1), 4] <- xTerm
        tot.mat[which(tot.mat[, 2] == maxDepth + 1), 5] <- y.heights[maxDepth + 1]
      }
      if (shape1 == "rad") {
        tot.mat[which(tot.mat[, 2] == maxDepth + 1), 4] <- theta
        wch <- which(tot.mat[, 2] == maxDepth + 1)
        tot.mat[wch[seq_along(wch) %% alt == 0], 5] <- r[i]
        tot.mat[wch[seq_along(wch) %% alt == 1], 5] <- r[i + 1]
        tot.mat[wch[seq_along(wch) %% alt == 2], 5] <- r[i + 2]
      }
      for (s in unique(tot.mat[which(tot.mat[, 2] == i), 3])) {
        cond1 <- which((tot.mat[, 2] == (i - 1)) & (tot.mat[, 3] == s))
        cond2 <- which((tot.mat[, 3] - s) == min(tot.mat[cond1, 3] - s))
        final.cond <- cond2[cond2 %in% cond1]
        pointer <- tot.mat[final.cond, 1]
        tot.mat[which((tot.mat[, 2] == i) & (tot.mat[, 3] == s)), 6] <- pointer
      }
    } else {
      for (s in unique(tot.mat[which(tot.mat[, 2] == i), 3])) {
        if (i != 1) {
          cond1 <- which((tot.mat[, 2] == (i - 1)) & (tot.mat[, 3] <= s))
          cond1 <- max(cond1)
          if (0 %in% (tot.mat[cond1, 3] - s)) {
            cond2 <- which((tot.mat[, 3] - s) == 0)
          } else {
            cond2 <- which((tot.mat[, 3] - s) == min(tot.mat[cond1, 3] - s))
          }
          final.cond <- cond2[cond2 %in% cond1]
          pointer <- tot.mat[final.cond, 1]
          tot.mat[which((tot.mat[, 2] == i) & (tot.mat[, 3] == s)), 6] <- pointer
        }
        if (i == 1) {
          tot.mat[which(tot.mat[, 2] == 1), 6] <- 0
        }

        if (shape1 == "trad") {
          c.row <- tot.mat[which((tot.mat[, 2] == i) & (tot.mat[, 3] == s)), ] # current row being editted
          x.val <- mean(tot.mat[which(tot.mat[, 6] == c.row[1]), 4])
          y.val <- y.heights[i]
          tot.mat[which((tot.mat[, 2] == i) & (tot.mat[, 3] == s)), 4] <- x.val
          tot.mat[which((tot.mat[, 2] == i) & (tot.mat[, 3] == s)), 5] <- y.val
        }
        if (shape1 == "rad") {
          c.row <- tot.mat[which((tot.mat[, 2] == i) & (tot.mat[, 3] == s)), ] # current row being editted
          theta.val <- mean(tot.mat[which(tot.mat[, 6] == c.row[1]), 4])
          r.val <- r[i]
          tot.mat[which((tot.mat[, 2] == i) & (tot.mat[, 3] == s)), 4] <- theta.val
          tot.mat[which((tot.mat[, 2] == i) & (tot.mat[, 3] == s)), 5] <- r.val
        }
      }
    }
  }

  ##################### Terminal node height alternator ########################
  if (alt > 1 && shape1 == "trad") {
    indeces.alt <- which(tot.mat[, 2] == maxDepth + 1)
    if (alt == 2) {
      if (length(indeces.alt) %% 2 == 0) {
        tot.mat[indeces.alt, 5] <- c(y.heights[maxDepth + 1], y.heights[maxDepth + 2])
      } else {
        tot.mat[indeces.alt, 5] <- c(
          rep(c(
            y.heights[maxDepth + 1],
            y.heights[maxDepth + 2]
          ), floor(length(indeces.alt) / 2)),
          y.heights[maxDepth + 1]
        )
      }
    }
    if (alt == 3) {
      if (length(indeces.alt) %% 3 == 0) {
        tot.mat[indeces.alt, 5] <- c(y.heights[maxDepth + 1], y.heights[maxDepth + 2], y.heights[maxDepth + 3])
      }
      if (length(indeces.alt) %% 3 == 1) {
        tot.mat[indeces.alt, 5] <- c(
          rep(c(
            y.heights[maxDepth + 1],
            y.heights[maxDepth + 2],
            y.heights[maxDepth + 3]
          ), floor(length(indeces.alt) / 3)),
          y.heights[maxDepth + 1]
        )
      }
      if (length(indeces.alt) %% 3 == 2) {
        tot.mat[indeces.alt, 5] <- c(
          rep(c(
            y.heights[maxDepth + 1],
            y.heights[maxDepth + 2],
            y.heights[maxDepth + 3]
          ), floor(length(indeces.alt) / 3)),
          y.heights[maxDepth + 1], y.heights[maxDepth + 2]
        )
      }
    }
  }

  ############################## Plotting lines ################################
  plot(NA, NA, xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), axes = FALSE)
  if (shape1 == "trad") {
    for (k in dim(tot.mat)[1]:1) {
      if (tot.mat[k, 6] > 1) {
        if (shape2 == "right") {
          segments(tot.mat[k, 4], tot.mat[k, 5],
            tot.mat[k, 4],
            tot.mat[which(tot.mat[, 1] == tot.mat[k, 6]), 5],
            col = if (k < maxDepth + 2) "black" else "black", lwd = lwd
          ) # vertical line
          segments(tot.mat[k, 4], tot.mat[which(tot.mat[, 1] == tot.mat[k, 6]), 5],
            tot.mat[which(tot.mat[, 1] == tot.mat[k, 6]), 4],
            tot.mat[which(tot.mat[, 1] == tot.mat[k, 6]), 5],
            col = if (k < maxDepth + 2) "black" else "black", lwd = lwd
          ) # horizontal line
        } else {
          segments(tot.mat[k, 4], tot.mat[k, 5],
            tot.mat[which(tot.mat[, 1] == tot.mat[k, 6]), 4],
            tot.mat[which(tot.mat[, 1] == tot.mat[k, 6]), 5],
            col = if (k < maxDepth + 2) "black" else "black", lwd = lwd
          )
        }
      }
    }
  }
  if (shape1 == "rad") {
    for (k in dim(tot.mat)[1]:1) {
      if (tot.mat[k, 6] > 1) {
        new.theta <- tot.mat[which(tot.mat[, 1] == tot.mat[k, 6]), 4]
        new.r <- tot.mat[which(tot.mat[, 1] == tot.mat[k, 6]), 5]
        segments(tot.mat[k, 5] * cos(tot.mat[k, 4]), tot.mat[k, 5] * sin(tot.mat[k, 4]),
          new.r * cos(tot.mat[k, 4]), new.r * sin(tot.mat[k, 4]),
          col = if (k < maxDepth + 2) "black" else "black", lwd = lwd
        ) # radial line
        lines(new.r * cos(seq(tot.mat[k, 4], new.theta, len = 100)),
          new.r * sin(seq(tot.mat[k, 4], new.theta, len = 100)),
          col = if (k < maxDepth + 2) "black" else "black", lwd = lwd
        )
      }
    }
  }

  ########################## Plotting words and boxes ##########################
  mywords <- c(na.omit(as.vector(t(var.mat.na))))
  mywords <- gsub("FILLER", "_", mywords)
  mywords <- gsub("FILLERPERIOD", "\\.", mywords)

  if (shape1 == "trad") {
    rect(tot.mat[-1, 4] - strwidth(mywords, "user") * cex / 1.9,
      tot.mat[-1, 5] - strheight(mywords, "user") * cex / 1.1,
      tot.mat[-1, 4] + strwidth(mywords, "user") * cex / 1.9,
      tot.mat[-1, 5] + strheight(mywords, "user") * cex / 1.1,
      col = col, border = col, lwd = lwd
    )
    lx <- tot.mat[-1, 4]
    ly <- tot.mat[-1, 5]
  }
  if (shape1 == "rad") {
    rect(tot.mat[-1, 5] * cos(tot.mat[-1, 4]) - strwidth(mywords, "user") * cex / 1.9,
      tot.mat[-1, 5] * sin(tot.mat[-1, 4]) - strheight(mywords, "user") * cex / 1.1,
      tot.mat[-1, 5] * cos(tot.mat[-1, 4]) + strwidth(mywords, "user") * cex / 1.9,
      tot.mat[-1, 5] * sin(tot.mat[-1, 4]) + strheight(mywords, "user") * cex / 1.1,
      col = col, border = col, lwd = lwd
    )
    lx <- tot.mat[-1, 5] * cos(tot.mat[-1, 4])
    ly <- tot.mat[-1, 5] * sin(tot.mat[-1, 4])
  }
  text(lx, ly, mywords, cex = cex, col = "black")
}
