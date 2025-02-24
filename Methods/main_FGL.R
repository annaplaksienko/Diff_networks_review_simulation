library(parallel)
library(jewel)
library(JGL)
library(dplyr)
library(igraph)

FGL_parallel <- function(data, lambda1, lambda2, thresh, 
                         true_diff_graph) {
    
    nlambda1 <- length(lambda1)
    nlambda2 <- length(lambda2)
    nthresh <- length(thresh)
    
    perf <- as.data.frame(matrix(0, nrow = nlambda1 * nlambda2 * nthresh, 
                                 ncol = 10))
    colnames(perf) <- c("lambda1", "lambda2", "threshold",
                        "TP", "TN", "FP", "FN", 
                        "FDR", "power",
                        "size")
    perf[, "lambda1"] <- rep(lambda1, each = nlambda2 * nthresh)
    perf[, "lambda2"] <- rep(rep(lambda2, each = nthresh), nlambda1)
    perf[, "threshold"] <- rep(thresh, nlambda1 * nlambda2)
    
    ind <- 1
    for (t in 1:nlambda1) {
        print(t)
        for (s in 1:nlambda2) {
            print(s)
            res_FGL <- JGL(data, penalty = "fused", 
                           lambda1 = lambda1[t], lambda2 = lambda2[s],
                           return.whole.theta = TRUE)$theta
            
            for (r in 1:nthresh) {
                curr_res <- abs(res_FGL[[1]] - res_FGL[[2]]) > thresh[r]
                
                diag(curr_res) <- 0
                perf[ind, "size"] <- sum(curr_res) / 2
                curr_perf <- suppressMessages(evaluatePerformance(true_diff_graph, 
                                                                  curr_res))
                perf[ind, c("TP", "TN", "FP", "FN")] <- curr_perf
                perf[ind, "FDR"] <- curr_perf["FP"] / 
                    pmax(curr_perf["TP"] + curr_perf["FP"], 1)
                perf[ind, "power"] <- curr_perf["TP"] / (sum(true_diff_graph) / 2) 
                
                ind <- ind + 1
            }    
        }
    }
    
    return(perf)
}

#we do not provide this file as it is too heavy. 
#Please produce it yourself using Generate_X.R
load("G_diff_X_all.Rdata")
#where to save the results
filename <- "results_FGL.Rdata"

ncores <- 8

lambda1 <- c(0.1, 0.25)
nlambda1 <- length(lambda1)
lambda2 <- c(seq(0.01, 0.1, length.out = 10),
             seq(0.1, 0.4, length.out = 11))
nlambda2 <- length(lambda2)

thresh <- c(10^-1, 10^-2, 10^-3, 10^-4, 10^-6)
nthresh <- length(thresh)

message("Running fused graphical lasso (FGL)")

perf_list <- perf_summarized_list <- vector(mode = "list", 
                                            length = length(G_diff_list) * 2)
names(perf_list) <- names(perf_summarized_list) <- paste(rep(names(G_diff_list), 
                                                             each = 2), 
                                                         rep(names(X_list[[1]]), 
                                                             length(names(G_diff_list))),
                                                         sep = "_")
ind <- 1

for (i in 1:length(G_diff_list)) {
    G_diff <- as.matrix(as_adjacency_matrix(G_diff_list[[i]]))
    timestamp(prefix = "@@-", suffix = "-@@")
    message(names(G_diff_list)[i])
    
    for (j in 1:length(X_list[[i]])) {
        X_reps <- X_list[[i]][[j]]
        nreps <- length(X_reps)
        timestamp(prefix = "*-", suffix = "-*")
        message(names(X_list[[i]])[j])
        
        cl <- makeCluster(ncores)
        clusterEvalQ(cl, library("JGL"))
        clusterEvalQ(cl, library("jewel"))
        results <- vector(mode = "list", length = nreps)
        names(results) <- paste("Realization", 1:nreps, sep = "_")
        results <- clusterApply(cl, X_reps, FGL_parallel, 
                                lambda1 = lambda1, 
                                lambda2 = lambda2, 
                                thresh = thresh,
                                true_diff_graph = G_diff)
        stopCluster(cl)
        
        perf <- do.call(rbind, results)
        perf_list[[ind]] <- data.frame("Realisation" = rep(1:nreps, 
                                                           each = nlambda1 * nlambda2 * nthresh), 
                           perf)
        
        perf_summarized_list[[ind]] <- perf %>%
            group_by(lambda1, lambda2, threshold) %>%
            summarize(across(TP:size, mean))
      
        ind <- ind + 1
        
        save(list = c("perf_list", "perf_summarized_list"), 
             file = filename)
    }
}

{
    settings <- as.data.frame(settings[rep(seq_len(nrow(settings)), each = 2), ])
    settings$sample_size <- rep(c("100 samples", "400 samples"),
                                dim(settings)[1] / 2)
    
    npar <- dim(perf_summarized_list[[1]])[1]
    perf <-  do.call(rbind, perf_summarized_list)
    perf$graph_type <- rep(settings$graph_type, each = npar)
    perf$true_diff_size <- rep(settings$true_diff_size, each = npar)
    perf$true_diff_size <- factor(perf$true_diff_size,
                                  levels = c("50", "100"))
    perf$true_G1_size <- rep(settings$true_G1_size, each = npar)
    perf$true_G1_size <- factor(perf$true_G1_size,
                                levels = c("200", "400"))
    perf$sample_size <- rep(settings$sample_size, each = npar)
    perf$method <- rep("FGL", dim(perf)[1])
    perf_FGL_full <- perf
    remove(perf)
    
    perf_FGL <- perf_FGL_full[perf_FGL_full$threshold == 0.001, ]
    perf_FGL <- perf_FGL[perf_FGL$lambda1 == 0.1, ]
    perf_FGL$param <- perf_FGL$lambda2
    perf_FGL$lambda1 <- perf_FGL$lambda2 <- perf_FGL$threshold <- NULL
    
    save(list = c("perf_list", "perf_summarized_list", 
                  "perf_FGL_full", "perf_FGL"), 
         file = filename)
}

message("Finished!")
timestamp(prefix = "#-", suffix = "-#")





