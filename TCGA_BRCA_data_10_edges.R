library(DiffGraph)
library(JGL)
library(DiffNetFDR)
library(DCA)
library(Grace)

#data from DiffNetFDR package
load("Real_data/TCGA.BRCA.rda")
X <- TCGA.BRCA$X[c(1, 3)]
names(X) <- c("basal", "lumA")

#you can either run the code for each method or load the saved results
#load(Real_data/Real_data_for_paper_10_edges.Rdata)

#JGL fused (FGL)
{
  res <- JGL(X, penalty = "fused", lambda1 = 0.1, lambda2 = 0.392,
             return.whole.theta = TRUE)$theta
  g_diff_est <- abs(res[[1]] - res[[2]]) > 10^(-3)
  diag(g_diff_est) <- 0
  res_FGL_10 <- graph_from_adjacency_matrix(g_diff_est, mode = "undirected")
  gsize(res_FGL_10)
  res_FGL_10 <- delete.vertices(res_FGL_10, which(degree(res_FGL_10) == 0))
}

#Dtrace
{
  res_Dtrace_10 <- Dtrace(X, lambda = 0.555)$Delta.graph.full
  gsize(res_Dtrace_10)
  
  res_Dtrace_10 <- delete.vertices(res_Dtrace_10, 
                                   which(degree(res_Dtrace_10) == 0))
}

#DiffNetFDR (difference in precision matrices)
{
  X_long <- do.call(rbind, X)
  labels <- c(rep(names(X)[1], dim(X[[1]])[1]), 
              rep(names(X)[2], dim(X[[2]])[1]))
  
  res_DiffNetFDR_10_precmat <- DiffNet.FDR(X_long, labels, alpha = 0.025, 
                                test.type = "pmat",
                                parallel = TRUE)$Diff.net.connected
  gsize(res_DiffNetFDR_10_precmat)
}

#DiffNetFDR (difference in partial correlations)
{
  res_DiffNetFDR_11_partcor <- DiffNet.FDR(X_long, labels, alpha = 0.06, 
                                   test.type = "pcor",
                                   parallel = TRUE)$Diff.net.connected
  gsize(res_DiffNetFDR_11_partcor)
}

#DCA
source("Methods/DCA_func_corrected.R")
{
  res_DCA_10 <- DCA_func(X, alpha = c(0.02))
  res_DCA_10 <- graph_from_adjacency_matrix(res_DCA_10[[1]], 
                                            mode = "undirected")
  gsize(res_DCA_10)
  
  res_DCA_10 <- delete_vertices(res_DCA_10, which(degree(res_DCA_10) == 0))
}

#plot union_10
{
  vsize <- 1
  vcol <- "white"
  vfrcol <- "white"
  vlcex <- 1.1
  ewidth <- 2.5
  par(mar = c(1, 1, 1, 1) + 0.5)
  
  union_10_multiple <- rbind(as_edgelist(res_FGL_10), 
                             as_edgelist(res_DiffNetFDR_10_precmat), 
                             as_edgelist(res_DiffNetFDR_11_partcor),
                             as_edgelist(res_DCA_10),
                             as_edgelist(res_Dtrace_10))
  
  union_10_multiple <- graph_from_edgelist(union_10_multiple, 
                                           directed = FALSE)
  
  color <- c(rep("#F8766D", gsize(res_FGL_10)),
             rep("#00BFC4", gsize(res_DiffNetFDR_10_precmat)),
             rep("#619CFF", gsize(res_DiffNetFDR_11_partcor)),
             rep("#F564E3", gsize(res_DCA_10)),
             rep("#00BA38", gsize(res_Dtrace_10)))
  
  
  #coords <- tkplot(union_10_multiple, 
                   #vertex.size = vsize, vertex.color = vcol, 
                   #vertex.frame.color = vfrcol, vertex.label.cex = vlcex,
                   #edge.width = 3, edge.color = color)
  #layout <- tk_coords(coords)
  
  plot.igraph(union_10_multiple, layout = layout, 
              vertex.size = vsize, vertex.color = vcol, 
              vertex.frame.color = vfrcol, vertex.label.cex = vlcex,
              edge.width = 3, edge.color = color)

}

save(list = c("X", "res_FGL_10", "res_Dtrace_10", 
              "res_DiffNetFDR_10_precmat", 
              "res_DiffNetFDR_11_partcor",
              "res_DCA_10",
              "union_10_multiple",
              "layout", "color"), 
     file = "Real_data_for_paper_10_edges.Rdata")

