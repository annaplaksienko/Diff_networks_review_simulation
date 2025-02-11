library(parallel)
library(dplyr)
library(ONDSA)
library(igraph)
library(jewel)

ONSDA_parallel <- function(data, alpha_vec, true_diff_graph) {

    nalpha <- length(alpha_vec)
    perf <- as.data.frame(matrix(0, nrow = nalpha, ncol = 8))
    colnames(perf) <- c("alpha",
                        "TP", "TN", "FP", "FN", 
                        "FDR", "power",
                        "size")
    perf[ , "alpha"] <- alpha_vec
    
    # Standardize the raw data for each group
    Data_1Std <- apply(data$X1, 2, standardize)
    Data_2Std <- apply(data$X2, 2, standardize)
    
    # Set parameters
    p <- dim(Data_1Std)[2]
    K <- 2
    N <- c(dim(Data_1Std)[1], dim(Data_2Std)[1])
    colnames(Data_1Std) <- colnames(Data_2Std) <- varnames <- V(true_diff_graph)$name <-
        paste("X", 1:p, sep = "")
    
    true_diff_graph_matrix <- as.matrix(as_adjacency_matrix(true_diff_graph))
    
    
    # Estimate group-specific precision matrices using FastGGM
    # lambda_value is an array corresponding to lambda of each group
    lambda_value <- sqrt(2 * log(p / sqrt(N)) / N)
    fastggm_1 <- FastGGM::FastGGM_Parallel(Data_1Std, lambda = lambda_value[1])
    fastggm_2 <- FastGGM::FastGGM_Parallel(Data_2Std, lambda = lambda_value[2])
    
    # Create list of precision matrices
    Omega <- list(fastggm_1$precision, fastggm_2$precision)


    for (t in 1:nalpha) {
        # Run ONDSA with the estimated precision matrices
        # modification could have been made for ONDSA to take as an input 
        # a vector of alpha values but the method is fast enough that we skipped this
        curr_alpha <- alpha_vec[t]
        result <- ONDSA(Omega, p, K, N, curr_alpha, varnames)
        
        est_diff_graph <- graph_from_edgelist(as.matrix(result$differential_structures[, 3:4]),
                                              directed = FALSE)
        perf[t, "size"] <- gsize(est_diff_graph)
        
        est_diff_graph_matrix <- as.matrix(as_adjacency_matrix(est_diff_graph))
        p_est <- dim(est_diff_graph_matrix)[1]
        
        if (p_est == 0) {
            est_diff_graph_matrix_full <- matrix(0, nrow = p, ncol = p)
            colnames(est_diff_graph_matrix_full) <- rownames(est_diff_graph_matrix_full) <- 
                colnames(true_diff_graph_matrix)
        } else if (p_est != p) {
            isolated_vertices <- colnames(true_diff_graph_matrix)[!colnames(true_diff_graph_matrix) %in%
                                                                      colnames(est_diff_graph_matrix)]
            
            est_diff_graph_matrix_full <- matrix(0, nrow = p, ncol = p)
            est_diff_graph_matrix_full[1:p_est, 1:p_est] <- est_diff_graph_matrix
            colnames(est_diff_graph_matrix_full) <- rownames(est_diff_graph_matrix_full) <- 
                c(colnames(est_diff_graph_matrix), isolated_vertices)
            
            est_diff_graph_matrix_full <- est_diff_graph_matrix_full[colnames(true_diff_graph_matrix), 
                                                                     colnames(true_diff_graph_matrix)]
        } else {
            est_diff_graph_matrix_full <- est_diff_graph_matrix
            est_diff_graph_matrix_full <- est_diff_graph_matrix_full[colnames(true_diff_graph_matrix), 
                                                                     colnames(true_diff_graph_matrix)]
        }
        
        #curr_perf <- evaluatePerformance_graphs(true_diff_graph, est_diff_graph)
        
        #perf[t, c("TP", "TN", "FP", "FN", "FDR", "power")] <- curr_perf
        #perf[t, "size"] <- gsize(est_diff_graph)
        
        curr_perf <- suppressMessages(evaluatePerformance(true_diff_graph_matrix, 
                                                          est_diff_graph_matrix_full))
        perf[t, c("TP", "TN", "FP", "FN")] <- curr_perf
        perf[t, "FDR"] <- curr_perf["FP"] / 
            pmax(curr_perf["TP"] + curr_perf["FP"], 1)
        perf[t, "power"] <- curr_perf["TP"] / (sum(true_diff_graph_matrix) / 2)
    }
    
    return(perf)
}

#we do not provide this file as it is too heavy. Please produce it yourself usign Generate_X.R
load("G_X_full.Rdata")

ncores <- 50
alpha <-  c(0.01, 0.02, 0.03, 0.04,
            seq(0.05, 0.995, by = 0.025))
nalpha <- length(alpha)

message("Running ONSDA")

perf_list <- perf_summarized_list <- vector(mode = "list", 
                                            length = length(G_diff_list) * 2)
names(perf_list) <- names(perf_summarized_list) <- paste(rep(names(G_diff_list), 
                                                             each = 2), 
                                                         rep(names(X_list[[1]]), 
                                                             length(names(G_diff_list))),
                                                         sep = "_")
ind <- 1

for (i in 1:length(G_diff_list)) {
    G_diff <- G_diff_list[[i]]
    timestamp(prefix = "@@-", suffix = "-@@")
    message(names(G_diff_list)[i])
    
    for (j in 1:length(X_list[[i]])) {
        X_reps <- X_list[[i]][[j]]
        nreps <- length(X_reps)
        timestamp(prefix = "*-", suffix = "-*")
        message(names(X_list[[i]])[j])
        
        cl <- makeCluster(ncores)
        clusterEvalQ(cl, library("ONDSA"))
        clusterEvalQ(cl, library("FastGGM"))
        clusterEvalQ(cl, library("igraph"))
        clusterEvalQ(cl, library("jewel"))
        results <- vector(mode = "list", length = nreps)
        names(results) <- paste("Realization", 1:nreps, sep = "_")
        results <- clusterApply(cl, X_reps, ONSDA_parallel, 
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
             file = "results_ONDSA.Rdata")
    }
}

message("Finished!")
timestamp(prefix = "#-", suffix = "-#")

#construct performance table
{
    settings <- settings[rep(seq_len(nrow(settings)), each = 2), ]
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
    perf$method <- rep("ONDSA", dim(perf)[1])
    perf$param <- perf$alpha
    perf$alpha <- NULL
    perf_ONDSA <- perf
    
    save(list = c("perf_list", "perf_summarized_list", "perf_ONDSA"), 
         file = "results_ONDSA.Rdata")
}