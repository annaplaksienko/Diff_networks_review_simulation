setwd("/Users/annapla/Documents/2_OCBE/3_Code/Diff_networks_review_simulation/")

library(DiffGraph)
library(JGL)
library(DiffNetFDR)
library(DCA)
library(Grace)
library(matrixcalc)
library(ONDSA)

#data from DiffNetFDR package
load("Real_data/TCGA.BRCA.rda")
X <- TCGA.BRCA$X[c(1, 3)]
names(X) <- c("basal", "lumA")
p <- dim(X$basal)[2]

#you can either run the code for each method or load the saved results
#load(Real_data/Real_data_for_paper_10_edges.Rdata)

#JGL fused (FGL) AIC
#produces much denser network than all other methods so hard to compare
{
    lambda1_grid <- seq(0.001, 1, length.out = 50)
    lambda2_grid <- seq(0.01, 0.9, length.out = 30)
    
    AIC_vec_step1 <- rep(NA, length(lambda1_grid))
    AIC_vec_step2 <- rep(NA, length(lambda2_grid))

    n1 <- dim(X$basal)[1]
    n2 <- dim(X$lumA)[1]
    Sigma1 <- cov(X$basal)
    Sigma2 <- cov(X$lumA)  
    
    for (l in 1:length(lambda1_grid)) {
        res <- JGL(X, penalty = "fused", lambda1 = lambda1_grid[l], 
                   lambda2 = lambda2_grid[1],
                   return.whole.theta = TRUE)$theta
        Theta1 <- Theta1_zero_diag <- res[[1]]
        Theta2 <- Theta2_zero_diag <- res[[2]]
        
        diag(Theta1_zero_diag) <- diag(Theta2_zero_diag) <- 0
        
        AIC_vec_step1[l] <- n1 * (sum(diag(Sigma1 %*% Theta1)) - 
                                      log(det(Theta1))) + 
            2 * sum(Theta1_zero_diag != 0) +
            n2 * (sum(diag(Sigma2 %*% Theta2)) - 
                      log(det(Theta2))) + 
            2 * sum(Theta2_zero_diag != 0)
    }
    
    lambda1_opt <- lambda1_grid[which.min(AIC_vec_step1)]
    plot(AIC_vec_step1, main = "AIC step 1")
    
    for (l in 1:length(lambda2_grid)) {
        res <- JGL(X, penalty = "fused", lambda1 = lambda1_opt, 
                   lambda2 = lambda2_grid[l],
                   return.whole.theta = TRUE)$theta
        Theta1 <- Theta1_zero_diag <- res[[1]]
        Theta2 <- Theta2_zero_diag <- res[[2]]
        
        diag(Theta1_zero_diag) <- diag(Theta2_zero_diag) <- 0
        
        AIC_vec_step2[l] <- n1 * (sum(diag(Sigma1 %*% Theta1)) - 
                                      log(det(Theta1))) + 
            2 * sum(Theta1_zero_diag != 0) +
            n2 * (sum(diag(Sigma2 %*% Theta2)) - 
                      log(det(Theta2))) + 
            2 * sum(Theta2_zero_diag != 0)
    }
    
    plot(AIC_vec_step2, main = "AIC step 2")
    
    lambda2_opt <- lambda2_grid[which.min(BIC_vec_step2)]
    
    #as minimum is achieved at the end of the interval, it is not the real minimum.
    #we suggest to consider the "elbow", i.e. where curve slightly changes its shape
    
    opt_res <- JGL(X, penalty = "fused", lambda1 = lambda1_opt, 
                            lambda2 = lambda2_grid[10],
                            return.whole.theta = TRUE)$theta
    g_diff_est <- abs(opt_res[[1]] - opt_res[[2]]) > 10^(-3)
    diag(g_diff_est) <- 0
    res_FGL_opt <- graph_from_adjacency_matrix(g_diff_est, mode = "undirected")
    gsize(res_FGL_opt)
    
}

#JGL fused (FGL) BIC
#produces a more sparse network than the one with AIC
#note, however, that both for AIC and BIC we didn't find a true minimum for lambda2
{
    lambda1_grid <- seq(0.001, 1, length.out = 50)
    lambda2_grid <- seq(0.01, 0.9, length.out = 30)

    BIC_vec_step1 <- rep(NA, length(lambda1_grid))
    BIC_vec_step2 <- rep(NA, length(lambda2_grid))
    
    n1 <- dim(X$basal)[1]
    n2 <- dim(X$lumA)[1]
    Sigma1 <- cov(X$basal)
    Sigma2 <- cov(X$lumA)  
    
    for (l in 1:length(lambda1_grid)) {
        res <- JGL(X, penalty = "fused", lambda1 = lambda1_grid[l], 
                   lambda2 = lambda2_grid[1],
                   return.whole.theta = TRUE)$theta
        Theta1 <- Theta1_zero_diag <- res[[1]]
        Theta2 <- Theta2_zero_diag <- res[[2]]
        
        diag(Theta1_zero_diag) <- diag(Theta2_zero_diag) <- 0
        
        BIC_vec_step1[l] <- n1 * (sum(diag(Sigma1 %*% Theta1)) - 
                                      log(det(Theta1))) + 
            log(n1) * sum(Theta1_zero_diag != 0) +
            n2 * (sum(diag(Sigma2 %*% Theta2)) - 
                      log(det(Theta2))) + 
            log(n2) * sum(Theta2_zero_diag != 0)
    }
    
    lambda1_opt_BIC <- lambda1_grid[which.min(BIC_vec_step1)]
    plot(BIC_vec_step1, main = "BIC step 1")
    
    for (l in 1:length(lambda2_grid)) {
        res <- JGL(X, penalty = "fused", lambda1 = lambda1_opt_BIC, 
                   lambda2 = lambda2_grid[l],
                   return.whole.theta = TRUE)$theta
        Theta1 <- Theta1_zero_diag <- res[[1]]
        Theta2 <- Theta2_zero_diag <- res[[2]]
        
        diag(Theta1_zero_diag) <- diag(Theta2_zero_diag) <- 0
        
        BIC_vec_step2[l] <- n1 * (sum(diag(Sigma1 %*% Theta1)) - 
                                      log(det(Theta1))) + 
            log(n1) * sum(Theta1_zero_diag != 0) +
            n2 * (sum(diag(Sigma2 %*% Theta2)) - 
                      log(det(Theta2))) + 
            log(n2) * sum(Theta2_zero_diag != 0)
    }

    plot(BIC_vec_step2, main = "BIC step 2")
    #as minimum is achieved at the end of the interval, it is not the real minimum.
    #we suggest to consider the "elbow", i.e. where curve slightly changes its shape
    #lambda2_opt_BIC <- lambda2_grid[10]
    
    BIC_res <- JGL(X, penalty = "fused", lambda1 = lambda1_opt_BIC, 
                            lambda2 = lambda2_opt_BIC,
                            return.whole.theta = TRUE)$theta
    g_diff_est <- abs(BIC_res[[1]] - BIC_res[[2]]) > 10^(-3)
    diag(g_diff_est) <- 0
    res_FGL_BIC <- graph_from_adjacency_matrix(g_diff_est, mode = "undirected")
    gsize(res_FGL_BIC)
    
}

#Dtrace
{
    lambda_grid <- c(seq(0.001, 0.1, by = 0.005),
                     seq(0.11, 0.8, by = 0.01))
    BIC_vec <- rep(NA, length(lambda_grid))
    
    n1 <- dim(X$basal)[1]
    n2 <- dim(X$lumA)[1]
    Sigma1 <- cov(X$basal)
    Sigma2 <- cov(X$lumA)
    
    for (i in 1:length(lambda_grid)) {
        res_Dtrace <- Dtrace(X, lambda = lambda_grid[i])
        #values of the difference
        Delta <- res_Dtrace$Delta
        #support of the difference
        Delta_graph <- res_Dtrace$Delta.graph.full
        
        BIC <- (n1 + n2) * frobenius.norm(0.5 * (Sigma1 %*% Delta %*% Sigma2 + 
                                        Sigma2 %*% Delta %*% Sigma1) -
                                        Sigma1 + Sigma2) +
            log(n1 + n2) * gsize(Delta_graph)
        BIC_vec[i] <- BIC
    }
  
lambda_opt <- lambda_grid[which.min(BIC_vec)] 
res_Dtrace <- Dtrace(X, lambda = lambda_opt)
gsize(res_Dtrace$Delta.graph.full)
plot(BIC_vec)  
plot(BIC_vec[50:90])  
    
res_Dtrace_BIC <- res_Dtrace$Delta.graph.connected
}

#DiffNetFDR (difference in precision matrices)
{
  X_long <- do.call(rbind, X)
  labels <- c(rep(names(X)[1], dim(X[[1]])[1]), 
              rep(names(X)[2], dim(X[[2]])[1]))
  
  res_DiffNetFDR_0.05_precmat <- DiffNet.FDR(X_long, labels, alpha = 0.05, 
                                test.type = "pmat",
                                parallel = TRUE)$Diff.net.connected
  gsize(res_DiffNetFDR_0.05_precmat)
  
  res_DiffNetFDR_0.01_precmat <- DiffNet.FDR(X_long, labels, alpha = 0.01, 
                                             test.type = "pmat",
                                             parallel = TRUE)$Diff.net.connected
  gsize(res_DiffNetFDR_0.01_precmat)
}

#DiffNetFDR (difference in partial correlations)
{
  res_DiffNetFDR_0.05_partcor <- DiffNet.FDR(X_long, labels, alpha = 0.05, 
                                   test.type = "pcor",
                                   parallel = TRUE)$Diff.net.connected
  gsize(res_DiffNetFDR_0.05_partcor)
  
  res_DiffNetFDR_0.01_partcor <- DiffNet.FDR(X_long, labels, alpha = 0.01, 
                                             test.type = "pcor",
                                             parallel = TRUE)$Diff.net.connected
  gsize(res_DiffNetFDR_0.01_partcor)
}

#DCA
source("Methods/DCA_func_corrected.R")
{
  res_DCA <- DCA_func(X, alpha = c(0.05, 0.01))
  res_DCA_0.05 <- graph_from_adjacency_matrix(res_DCA$`alpha = 0.05`, 
                                            mode = "undirected")
  gsize(res_DCA_0.05)
  res_DCA_0.01 <- graph_from_adjacency_matrix(res_DCA$`alpha = 0.01`, 
                                              mode = "undirected")
  gsize(res_DCA_0.01)
  
  
  res_DCA_10 <- delete_vertices(res_DCA_10, which(degree(res_DCA_10) == 0))
}

#ONDSA
{
    # Standardize the raw data for each group
    Data_1Std <- apply(X$basal, 2, standardize)
    Data_2Std <- apply(X$lumA, 2, standardize)
    
    # Set parameters
    p <- dim(Data_1Std)[2]
    K <- 2
    N <- c(dim(Data_1Std)[1], dim(Data_2Std)[1])
    colnames(Data_1Std) <- colnames(Data_2Std) <- varnames <- 
        colnames(X$basal)
    
    # Estimate group-specific precision matrices using FastGGM
    # lambda_value is an array corresponding to lambda of each group
    lambda_value <- sqrt(2 * log(p / sqrt(N)) / N)
    fastggm_1 <- FastGGM::FastGGM_Parallel(Data_1Std, lambda = lambda_value[1])
    fastggm_2 <- FastGGM::FastGGM_Parallel(Data_2Std, lambda = lambda_value[2])
    
    # Create list of precision matrices
    Omega <- list(fastggm_1$precision, fastggm_2$precision)
    
    # Run ONDSA with the estimated precision matrices
    result <- ONDSA(Omega, p, K, N, alpha = 0.05, varnames)
    res_ONDSA_0.05 <- graph_from_edgelist(as.matrix(result$differential_structures[, 3:4]),
                                        directed = FALSE)
    gsize(res_ONDSA_0.05)
    
    res_ONDSA_0.05 <- delete_vertices(res_ONDSA_0.05, 
                                    which(degree(res_ONDSA_0.05) == 0))
    
}

#plot union_opt
{
  vsize <- 1
  vcol <- "white"
  vfrcol <- "white"
  vlcex <- 1.1
  ewidth <- 2.5
  par(mar = c(1, 1, 1, 1) + 0.5)
  
  union_opt <- rbind(as_edgelist(res_FGL_BIC), 
                             as_edgelist(res_DiffNetFDR_0.05_precmat), 
                             as_edgelist(res_DiffNetFDR_0.05_partcor),
                             as_edgelist(res_DCA_0.05),
                             as_edgelist(res_Dtrace_BIC),
                             as_edgelist(res_ONDSA_0.05))
  
  union_opt <- graph_from_edgelist(union_opt, directed = FALSE)
  
  color <- c(rep("#F8766D", gsize(res_FGL_BIC)),
             rep("#00BFC4", gsize(res_DiffNetFDR_0.05_precmat)),
             rep("#619CFF", gsize(res_DiffNetFDR_0.05_partcor)),
             rep("#F564E3", gsize(res_DCA_0.05)),
             rep("#00BA38", gsize(res_Dtrace_BIC)),
             rep("#B79F00", gsize(res_ONDSA_0.05)))
  
  coords <- tkplot(union_opt, 
                   vertex.size = vsize, vertex.color = vcol, 
                   vertex.frame.color = vfrcol, vertex.label.cex = vlcex,
                   edge.width = 3, edge.color = color)
  #layout <- tk_coords(coords)
  
  plot.igraph(union_opt, layout = layout, 
              vertex.size = vsize, vertex.color = vcol, 
              vertex.frame.color = vfrcol, vertex.label.cex = vlcex,
              edge.width = 3, edge.color = color)
  
  #generating legend
  test <- data.frame(x = c(1:6), y = c(1:6), 
                     method = c("FGL", "Testing prec mat", 
                                'Testing part corr', 
                                "ONDSA", "DCA", "Dtrace"))
  test$method <- factor(test$method, levels = 
                            c("FGL", "Testing prec mat", 
                              'Testing part corr', 
                              "ONDSA", "DCA", "Dtrace"))
  
  ggplot(test, aes(x = x, y = y, color = method)) + 
      geom_line(linewidth = 1.1) +  
      theme(text = element_text(size = 28), 
            plot.title = element_text(size = 24)) +
      scale_color_manual(values = c("#F8766D", "#00BFC4", 
                                    "#619CFF", "#B79F00", "#F564E3",
                                    "#00BA38"))
}

#about the graph
{
    gorder(union_opt)
    gsize(simplify(union_opt))
    table(count_multiple(union_opt)) 
}

save(list = c("X", "res_Dtrace_BIC", "BIC_vec", "lambda_grid", 
              "res_FGL_BIC", "lambda1_opt_BIC", "lambda2_opt_BIC",
              "res_DiffNetFDR_0.05_precmat", "res_DiffNetFDR_0.01_precmat", 
              "res_DiffNetFDR_0.05_partcor", "res_DiffNetFDR_0.01_partcor",
              "res_DCA_0.01", "res_DCA_0.05",
              "res_ONDSA_0.05", 
              "union_opt", "layout"), 
     file = "Real_data_for_paper_optimal_parameter.Rdata")

