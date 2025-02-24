library(parallel)
library(jewel)
library(dplyr)
library(igraph)
library(DCA)

DCA_parallel <- function(data, alpha_vec, true_diff_graph) {
  alpha_res <- vector(mode = "list", length = length(alpha_vec))
  names(alpha_res) <- paste("alpha", alpha_vec, sep = "_") 
  
  nalpha <- length(alpha_vec)
  perf <- as.data.frame(matrix(0, nrow = nalpha, ncol = 8))
  colnames(perf) <- c("param",
                      "TP", "TN", "FP", "FN", 
                      "FDR", "power",
                      "size")
  perf[ , "param"] <- alpha_vec
  
  #res_DCA <- DCA_func(data, alpha = alpha_vec)$diffnet_list
  res_DCA <- DCA_func(data, alpha = alpha_vec)
  
  for (t in 1:nalpha) {
    curr_res <- res_DCA[[t]]
    perf[t, "size"] <- sum(curr_res) / 2
    curr_perf <- suppressMessages(evaluatePerformance(true_diff_graph, 
                                                      curr_res))
    perf[t, c("TP", "TN", "FP", "FN")] <- curr_perf
    perf[t, "FDR"] <- curr_perf["FP"] / 
      pmax(curr_perf["TP"] + curr_perf["FP"], 1)
    perf[t, "power"] <- curr_perf["TP"] / (sum(true_diff_graph) / 2)
  }
  
  return(perf)
}

source("DCA_func_corrected.R")

#we do not provide this file as it is too heavy. 
#Please produce it yourself using Generate_X.R
load("G_diff_X_all.Rdata")
#for saving results
filename <- "results_DCA.Rdata"

message("Running DCA")

ncores <- 9
alpha <- c(0.0001, 0.005, 0.05, 
           seq(0.1, 0.95, by = 0.05),
           0.97, 0.99)
nalpha <- length(alpha)

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
        clusterEvalQ(cl, library("jewel"))
        clusterEvalQ(cl, library("DCA"))
        clusterEvalQ(cl, source("/Users/annapla/Documents/2_OCBE/3_Manuscript_1/Diff_networks_review_simulation/Methods/DCA_func_corrected.R"))
        results <- vector(mode = "list", length = nreps)
        results <- clusterApply(cl, X_reps, DCA_parallel, 
                                alpha_vec = alpha, 
                                true_diff_graph = G_diff)
        names(results) <- paste("Realization", 1:nreps, sep = "_")

        perf <- do.call(rbind, results)
        perf_list[[ind]] <- data.frame("Realisation" = rep(1:nreps, each = nalpha), 
                                       perf)
        
        perf_summarized_list[[ind]] <- perf %>%
          group_by(param) %>%
          summarize(across(TP:size, mean))
      
        ind <- ind + 1
        
        save(list = c("perf_list", "perf_summarized_list"), 
             file = filename)
        
        remove(results)
    }
}

{
    settings <- as.data.frame(settings[rep(seq_len(nrow(settings)), 
                                           each = 2), ])
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
    perf$method <- rep("DCA", dim(perf)[1])
    perf$param <- perf$alpha
    perf$alpha <- NULL
    perf_DCA <- perf
    remove(perf)
    
    save(list = c("perf_list", "perf_summarized_list", "perf_DCA"), 
         file = filename)
}

message("Finished!")
timestamp(prefix = "#-", suffix = "-#")

