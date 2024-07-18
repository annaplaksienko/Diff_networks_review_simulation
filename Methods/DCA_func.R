#this is a modification of the function from DCA package https://github.com/sen-zhao/DCA 
#allowing us to produce results for several alpha levels in one run

DCA_func <- function(X, alpha, rule = "OR") {
  
  test.method <- "GraceI"
  test.level <- "edge"
  m <- length(X)
  
  est <- vector(mode = "list", length = m)
  for (i in 1:m) {
    data <- X[[i]]
    n <- nrow(data)
    p <- ncol(data)
    adjacency <- matrix(0, nrow = p, ncol = p)
    for (j in 1:p) {
      lammin <- cv.glmnet(data[, -j], data[, j])$lambda.min
      b <- glmnet(data[, -j], data[, j], lambda = lammin)$beta
      index <- which(b != 0)
      if (length(index) != 0) {
        index[index >= j] <- index[index >= j] + 1
      }
      adjacency[index, j] <- 1
    }
    est[[i]] <- adjacency
  }
  
  common <- est[[1]]
  for (i in 1:m) {
    common <- common * est[[i]]
  }
  
  nalpha <- length(alpha)
  p.thresh <- 1 - (1 - alpha)^(1/m)
  
  diffnet <- vector(mode = "list", length = nalpha)
  names(diffnet) <- paste("alpha", alpha, sep = "_")
  for (l in 1:nalpha) {
    diffnet[[l]] <- matrix(0, ncol = p, nrow = p)
  }
  pval_mat <-  matrix(1, ncol = p, nrow = p)
  
  for (i in 1:m) {
    data <- X[[i]]
    n <- nrow(data)
    p <- ncol(data)
    for (j in 1:p) {
      MJC <- which(common[j, -j] == 0)
      if (length(MJC) != 0) {
        Graceres <- graceI.test(data[, j], data[, -j], 
                                lambda.2 = exp(seq(from = -5, to = 10, 
                                                   length.out = 30)),
                                verbose = FALSE)
        p.Grace <- Graceres$pvalue[MJC]
        MJC[MJC >= j] <- MJC[MJC >= j] + 1
        for (l in 1:nalpha) {
          diffnet[[l]][j, MJC] <- diffnet[[l]][j, MJC] + 
            (p.Grace < p.thresh[l])
          pval_mat[j, MJC] <- pmin(pval_mat[j, MJC], p.Grace)
        }
      }
    }
  }  
  
  for (l in 1:nalpha) {
    diffnet[[l]] <- diffnet[[l]] + t(diffnet[[l]])
    if (rule == "AND") {
      diffnet[[l]][diffnet[[l]] < 2] <- 0
      diffnet[[l]][diffnet[[l]] == 2] <- 1
    }
    else if (rule == "OR") {
      diffnet[[l]][diffnet[[l]] >= 1] <- 1
    }
  }

  return(list(diffnet_list = diffnet, pval_mat = pval_mat))
}