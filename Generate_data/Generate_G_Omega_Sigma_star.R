library(igraph)
#for is.positive.definite
library(matrixcalc)

source("ggcorrplot_with_limits.R")

#number of vertices
p <- 200

vsize <- 2
par(mar = c(1, 1, 1, 1) + 0.5)

#G1 has 200 edges
{
    #re-run until we have 2 hubs of similar sizes, each around 50
    {
        test <- barabasi.game(p, power = 1.7, directed = FALSE)
        sort(degree(test))
    }
    G1_graph <- test
    G1 <- as.matrix(as_adjacency_matrix(G1_graph))
    
    #make Omega1
    {
        Omega1 <- G1
        size <- sum(G1) / 2
        a <- 0.6
        b <- 0.9
        samp_right <- runif(ceiling(size / 2), min = a, max = b)
        samp_left <- runif(ceiling(size / 2), min = -b, max = -a)
        samp <- sample(c(samp_left, samp_right), size)
        ind <- 1
        for (i in 2:p) {
            for (j in 1:(i - 1)) {
                if (G1[i, j] != 0) {
                    Omega1[i, j] <- Omega1[j, i] <- samp[ind]
                    ind <- ind + 1
                }
            }
        }
        is.positive.definite(Omega1)
        diag(Omega1) <- rowSums(abs(Omega1)) + 0.1
        is.positive.definite(Omega1)
        ggcorrplot2(Omega1, limits = range(Omega1))
        Sigma1 <- solve(Omega1)
    }
    
    #for 50-200 remove 1 hub
    {
        #remove the vertex with the highest degree! here it was number 1
        G_diff_graph <- graph_from_edgelist(ends(G1_graph, 
                                                 incident(G1_graph, 1)), 
                                            directed = FALSE)
        gorder(G_diff_graph)
        gsize(G_diff_graph)
        G_diff <- as.matrix(as_adjacency_matrix(G_diff_graph))
        
        G2_graph <- igraph::difference(G1_graph, G_diff_graph)
        G2 <- as.matrix(as_adjacency_matrix(G2_graph))
        gsize(G2_graph)
        
        G_plot <- G1 + G2
        diag(G_plot) <- 0 
        G_plot <-  graph_from_adjacency_matrix(G_plot, mode = "undirected", 
                                               weighted = TRUE)
        E(G_plot)$color = 3 - get.edge.attribute(G_plot, "weight")
        
        plot(G_plot, vertex.frame.color = "white", 
             layout = layout, 
             vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
             edge.width = 3)
        
        E(G2_graph)$color <- 1
        plot(G2_graph, vertex.frame.color = "white", 
             layout = layout, 
             vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
             edge.width = 3)
        
        E(G_diff_graph)$color <- 2
        plot(G_diff_graph, vertex.frame.color = "white", 
             layout = layout, 
             vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
             edge.width = 3)
        
        Omega2 <- Omega1 * G2
        diag(Omega2) <- diag(Omega1)
        ggcorrplot2(Omega2, limits = range(Omega1))
        is.positive.definite(Omega2)
        Sigma2 <- solve(Omega2)
        
        save(list = c("G1", "G1_graph", 
                      "G2", "G2_graph", 
                      "G_diff", "G_diff_graph",
                      "Omega1", "Sigma1", "Omega2", "Sigma2",
                      "layout"),
             file = "G_Omega_Sigma_star_50_200_.Rdata")
        
    }
    
    #for 100-200 remove 2 hubs
    {
        #remove the vertex with the highest degree! here it was numbers 1 and 2
        G_diff_graph <- graph_from_edgelist(rbind(ends(G1_graph, 
                                                       incident(G1_graph, 1)),
                                                  ends(G1_graph, 
                                                       incident(G1_graph, 2))), 
                                            directed = FALSE)
        G_diff_graph <- simplify(G_diff_graph)
        gorder(G_diff_graph)
        gsize(G_diff_graph)
        G_diff <- as.matrix(as_adjacency_matrix(G_diff_graph))
        
        G2_graph <- igraph::difference(G1_graph, G_diff_graph)
        G2 <- as.matrix(as_adjacency_matrix(G2_graph))
        gsize(G2_graph)
        
        G_plot <- G1 + G2
        diag(G_plot) <- 0 
        G_plot <-  graph_from_adjacency_matrix(G_plot, mode = "undirected", 
                                               weighted = TRUE)
        E(G_plot)$color = 3 - get.edge.attribute(G_plot, "weight")
        
        plot(G_plot, vertex.frame.color = "white", 
             layout = layout, 
             vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
             edge.width = 3)
        
        E(G2_graph)$color <- 1
        plot(G2_graph, vertex.frame.color = "white", 
             layout = layout, 
             vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
             edge.width = 3)
        
        E(G_diff_graph)$color <- 2
        plot(G_diff_graph, vertex.frame.color = "white", 
             layout = layout, 
             vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
             edge.width = 3)
        
        Omega2 <- Omega1 * G2
        diag(Omega2) <- diag(Omega1)
        ggcorrplot2(Omega2, limits = range(Omega1))
        is.positive.definite(Omega2)
        Sigma2 <- solve(Omega2)
        ggcorrplot2(Sigma2, limits = range(Sigma2))
        
        save(list = c("G1", "G1_graph", 
                      "G2", "G2_graph", 
                      "G_diff", "G_diff_graph",
                      "Omega1", "Sigma1", "Omega2", "Sigma2",
                      "layout"),
             file = "G_Omega_Sigma_star_100_200_.Rdata")
        
    }
}

#G1 has 400 edges
{
    #re-run until we have similar hub sizes, each around 50
    {
        test <- barabasi.game(p, power = 1.5, m = 2, directed = FALSE)
        sort(degree(test))
        gsize(test)
    }
    G1_graph <- test
    G1 <- as.matrix(as_adjacency_matrix(G1_graph))
    
    #make Omega1
    {
        Omega1 <- G1
        size <- sum(G1) / 2
        a <- 0.6
        b <- 0.9
        samp_right <- runif(ceiling(size / 2), min = a, max = b)
        samp_left <- runif(ceiling(size / 2), min = -b, max = -a)
        samp <- sample(c(samp_left, samp_right), size)
        ind <- 1
        for (i in 2:p) {
            for (j in 1:(i - 1)) {
                if (G1[i, j] != 0) {
                    Omega1[i, j] <- Omega1[j, i] <- samp[ind]
                    ind <- ind + 1
                }
            }
        }
        is.positive.definite(Omega1)
        diag(Omega1) <- rowSums(abs(Omega1)) + 0.1
        is.positive.definite(Omega1)
        ggcorrplot2(Omega1)
        Sigma1 <- solve(Omega1)
        ggcorrplot2(Sigma1)
    }
    
    #for 50-400 remove 1 hub
    {
        #remove the vertex with the highest degree! here it was number 17
        G_diff_graph <- graph_from_edgelist(ends(G1_graph, 
                                                 incident(G1_graph, 17)), 
                                            directed = FALSE)
        gorder(G_diff_graph)
        gsize(G_diff_graph)
        G_diff <- as.matrix(as_adjacency_matrix(G_diff_graph))
        
        G2_graph <- igraph::difference(G1_graph, G_diff_graph)
        G2 <- as.matrix(as_adjacency_matrix(G2_graph))
        gsize(G2_graph)
        
        G_plot <- G1 + G2
        diag(G_plot) <- 0 
        G_plot <-  graph_from_adjacency_matrix(G_plot, mode = "undirected", 
                                               weighted = TRUE)
        E(G_plot)$color = 3 - get.edge.attribute(G_plot, "weight")
        
        #layout <- layout_nicely(G_plot)
        
        plot(G_plot, vertex.frame.color = "white", 
             layout = layout, 
             vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
             edge.width = 3)
        
        E(G2_graph)$color <- 1
        plot(G2_graph, vertex.frame.color = "white", 
             layout = layout, 
             vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
             edge.width = 3)
        
        E(G_diff_graph)$color <- 2
        plot(G_diff_graph, vertex.frame.color = "white", 
             layout = layout, 
             vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
             edge.width = 3)
        
        Omega2 <- Omega1 * G2
        diag(Omega2) <- diag(Omega1)
        ggcorrplot2(Omega2)
        is.positive.definite(Omega2)
        Sigma2 <- solve(Omega2)
        ggcorrplot2(Sigma2)
        
        save(list = c("G1", "G1_graph", 
                      "G2", "G2_graph", 
                      "G_diff", "G_diff_graph",
                      "Omega1", "Sigma1", "Omega2", "Sigma2",
                      "layout"),
             file = "G_Omega_Sigma_star_50_400_.Rdata")
        
    }
    
    #for 100-400 remove 2 hubs
    {
        #remove the vertices with the highest degree! here it was number 17 and 1
        G_diff_graph <- graph_from_edgelist(rbind(ends(G1_graph, 
                                                       incident(G1_graph, 17)),
                                                  ends(G1_graph, 
                                                       incident(G1_graph, 1))), 
                                            directed = FALSE)
        G_diff_graph <- simplify(G_diff_graph)
        gorder(G_diff_graph)
        gsize(G_diff_graph)
        G_diff <- as.matrix(as_adjacency_matrix(G_diff_graph))
        
        G2_graph <- igraph::difference(G1_graph, G_diff_graph)
        G2 <- as.matrix(as_adjacency_matrix(G2_graph))
        gsize(G2_graph)
        
        G_plot <- G1 + G2
        diag(G_plot) <- 0 
        G_plot <-  graph_from_adjacency_matrix(G_plot, mode = "undirected", 
                                               weighted = TRUE)
        E(G_plot)$color = 3 - get.edge.attribute(G_plot, "weight")
        
        plot(G_plot, vertex.frame.color = "white", 
             layout = layout, 
             vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
             edge.width = 3)
        
        E(G2_graph)$color <- 1
        plot(G2_graph, vertex.frame.color = "white", 
             layout = layout, 
             vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
             edge.width = 3)
        
        E(G_diff_graph)$color <- 2
        plot(G_diff_graph, vertex.frame.color = "white", 
             layout = layout, 
             vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
             edge.width = 3)
        
        Omega2 <- Omega1 * G2
        diag(Omega2) <- diag(Omega1)
        ggcorrplot2(Omega2)
        is.positive.definite(Omega2)
        Sigma2 <- solve(Omega2)
        ggcorrplot2(Sigma2)
        
        save(list = c("G1", "G1_graph", 
                      "G2", "G2_graph", 
                      "G_diff", "G_diff_graph",
                      "Omega1", "Sigma1", "Omega2", "Sigma2",
                      "layout"),
             file = "G_Omega_Sigma_star_100_400_.Rdata")
        
    }
}


