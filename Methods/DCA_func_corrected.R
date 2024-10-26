#this is a modification of the function from DCA package https://github.com/sen-zhao/DCA 
#correction for multiplicity of edge estimators
#it also allows us to produce results for several alpha levels in one run

diffcon <- function(X, alpha) {
  
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
    #using the OR rule
    diffnet[[l]][diffnet[[l]] >= 1] <- 1
  }
  
  return(list(diffnet_list = diffnet, pval_mat = pval_mat))
}

#adjustment of multiplicity for OR rule
adjustment <- function(diffnet, pvalues_mat, method = "BH") {
    pvalues_mat_symmetric <- 2 * pmin(pvalues_mat, t(pvalues_mat))
    
    res_list <- vector(mode = "list", length = length(diffnet))
    names(res_list) <- names(diffnet)
    for (l in 1:length(diffnet)) {
        diffnet_upp <- upper.tri(diffnet[[l]])
        temp <- which(diffnet_upp == 1, arr.ind = T)
        pval_v <- p.adjust(as.numeric(pvalues_mat_symmetric[temp]), 
                           method = method)
        res_p <- matrix(1, nrow = nrow(pvalues_mat), 
                        ncol = ncol(pvalues_mat))
        res_p[temp] <- pval_v
        res_list[[l]] <- pmin(res_p, t(res_p))
    }

    return(res_list)
}

diffnet_from_adjusted <- function(pval_mat_adj, alpha){
    
    diffnet_list <- vector(mode = "list", length = length(alpha))
    
    for (l in 1:length(alpha)) {
        diffnet <- pval_mat_adj[[l]] <= alpha[l]
        diffnet <- diffnet + t(diffnet)
        diffnet[diffnet >= 1] <- 1
        diffnet_list[[l]] <- diffnet
    }
    
    return(diffnet_list)
}

DCA_func <- function(X, alpha) {
    
    non_adj_res <- diffcon(X, alpha)
    
    adjusted_p_values <- adjustment(diffnet = non_adj_res$diffnet_list, 
                                    pvalues_mat = non_adj_res$pval_mat)
    
    final_res <- diffnet_from_adjusted(adjusted_p_values, alpha)
    names(final_res) <- paste("alpha = ", alpha, sep = "")
    
    return(final_res)
}

