library(parallel)
library(jewel)
library(DiffGraph)
library(dplyr)

DTrace_parallel <- function(data, lambda_vec, true_diff_graph) {
    lambda_res <- vector(mode = "list", length = length(lambda_vec))
    names(lambda_res) <- paste("lambda", lambda_vec, sep = "_") 
    
    nlambda <- length(lambda_vec)
    perf <- as.data.frame(matrix(0, nrow = nlambda, ncol = 8))
    colnames(perf) <- c("lambda",
                        "TP", "TN", "FP", "FN", 
                        "FDR", "power",
                        "size")
    perf[ , "lambda"] <- lambda_vec
    
    for (t in 1:nlambda) {
        res <- Dtrace(data, lambda = lambda_vec[t])$Delta.graph.full
        res <- as.matrix(as_adjacency_matrix(res))
        curr_perf <- suppressMessages(evaluatePerformance(true_diff_graph, 
                                                          res))
        
        perf[t, c("TP", "TN", "FP", "FN")] <- curr_perf
        perf[t, "FDR"] <- curr_perf["FP"] / 
            pmax(curr_perf["TP"] + curr_perf["FP"], 1)
        perf[t, "power"] <- curr_perf["TP"] / (sum(true_diff_graph) / 2)
        
        perf[t, "size"] <- sum(res) / 2
    }
    
    return(perf)
}

#we do not provide this file as it is too heavy. Please produce it yourself usign Generate_X.R
load("G_diff_X_all.Rdata")
#where to save the results
filename <- "results_DTrace.Rdata"

G_diff_list <- G_diff_list[-c(3), drop = FALSE]
X_list <- X_list[-c(3), drop = FALSE]
settings <- settings[-c(3), , drop = FALSE]

ncores <- 50
lambda <- seq(from = 0.08, to = 0.8, length.out = 40)
nlambda <- length(lambda)

message("Running DTrace, all settings")

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
        clusterEvalQ(cl, library("DiffGraph"))
        clusterEvalQ(cl, library("jewel"))
        results <- vector(mode = "list", length = nreps)
        names(results) <- paste("Realization", 1:nreps, sep = "_")
        results <- clusterApply(cl, X_reps, DTrace_parallel, 
                                lambda_vec = lambda, 
                                true_diff_graph = G_diff)
        stopCluster(cl)
        
        perf <- do.call(rbind, results)
        perf_list[[ind]] <- data.frame("Realisation" = rep(1:nreps, each = nlambda), 
                           perf)
        
        perf_summarized_list[[ind]] <- perf %>%
            group_by(lambda) %>%
            summarize(across(TP:size, mean))
      
        ind <- ind + 1
        
        save(list = c("perf_list", "perf_summarized_list"), 
             file = filename)
    }
}

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
perf$method <- rep("DTrace", dim(perf)[1])
perf$param <- perf$lambda
perf$lambda <- NULL
perf_DTrace <- perf
remove(perf)

save(list = c("perf_list", "perf_summarized_list", "perf_DTrace"), 
     file = filename)

message("Finished!")
timestamp(prefix = "#-", suffix = "-#")
