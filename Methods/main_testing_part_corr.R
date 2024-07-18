library(parallel)
library(jewel)
library(DiffNetFDR)
library(dplyr)

source("DiffNetFDR_functions.R")

Testing_part_corr_parallel <- function(data, alpha_vec, true_diff_graph) {
    alpha_res <- vector(mode = "list", length = length(alpha_vec))
    names(alpha_res) <- paste("alpha", alpha_vec, sep = "_") 
    
    nalpha <- length(alpha_vec)
    perf <- as.data.frame(matrix(0, nrow = nalpha, ncol = 8))
    colnames(perf) <- c("alpha",
                        "TP", "TN", "FP", "FN", 
                        "FDR", "power",
                        "size")
    perf[ , "alpha"] <- alpha_vec
    
    X_long <- do.call(rbind, data)
    labels <- c(rep(names(data)[1], dim(data[[1]])[1]), 
                rep(names(data)[2], dim(data[[2]])[1]))
    res_DiffNetFDR_precmat <- DiffNet.FDR_new(X_long, labels, 
                                              alpha = alpha_vec, 
                                              test.type = "pcor",
                                              parallel = FALSE)
    for (t in 1:nalpha) {
        curr_res <- res_DiffNetFDR_precmat[[t]]$Diff.net
        perf[t, "size"] <- gsize(curr_res)
        curr_perf <- suppressMessages(evaluatePerformance(true_diff_graph, 
                                                          as_adjacency_matrix(curr_res)))
        perf[t, c("TP", "TN", "FP", "FN")] <- curr_perf
        perf[t, "FDR"] <- curr_perf["FP"] / 
            pmax(curr_perf["TP"] + curr_perf["FP"], 1)
        perf[t, "power"] <- curr_perf["TP"] / (sum(true_diff_graph) / 2)
    }
    
    return(perf)
}

#we do not provide this file as it is too heavy. Please produce it yourself usign Generate_X.R
load("G_X_full.Rdata")

ncores <- 50
alpha <- c(0.001, 0.0025, 0.005, 
           0.01, 0.02, 0.03, 0.04,
           seq(0.05, 0.9, by = 0.025),
           seq(0.91, 0.995, by = 0.005))
nalpha <- length(alpha)

message("Running DiffNetFDR, testing for partial correlations equality")

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
        clusterEvalQ(cl, library("DiffNetFDR"))
        clusterEvalQ(cl, source("DiffNetFDR_functions.R"))
        clusterEvalQ(cl, library("jewel"))
        results <- vector(mode = "list", length = nreps)
        names(results) <- paste("Realization", 1:nreps, sep = "_")
        results <- clusterApply(cl, X_reps, Testing_part_corr_parallel, 
                                alpha_vec = alpha, 
                                true_diff_graph = G_diff)
        stopCluster(cl)
        
        perf <- do.call(rbind, results)
        perf_list[[ind]] <- data.frame("Realisation" = rep(1:nreps, each = nalpha), 
                           perf)
        
        perf_summarized_list[[ind]] <- perf %>%
            group_by(alpha) %>%
            summarize(across(TP:size, mean))
      
        ind <- ind + 1
        
        save(list = c("perf_list", "perf_summarized_list"), 
             file = "Results_testing_part_corr.Rdata")
    }
}

message("Finished!")
timestamp(prefix = "#-", suffix = "-#")



