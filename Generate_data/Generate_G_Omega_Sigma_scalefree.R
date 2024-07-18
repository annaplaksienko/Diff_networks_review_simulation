library(igraph)
#for is.positive.definite
library(matrixcalc)

source("ggcorrplot_with_limits.R")
#load the setting you need
#load("G1_Omega1_Sigma1_m2_400_power1_scalefree_random.Rdata")

vsize <- 2
par(mar = c(1, 1, 1, 1) + 0.5)

plot(G1_graph, vertex.frame.color = "white", 
     layout = layout, 
     vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
     edge.width = 3)

p <-  200
L <-  4
p_sub <-  p / L

#for 50-200 and 100-400 we remove one block 
{
    G2 <- G1
    G2[151:200, 151:200] <- 0
    ggcorrplot2(G2)
}

#for 100-200 we remove two blocks
{
    G2 <- G1
    G2[101:200, 101:200] <- 0
    ggcorrplot2(G2)
}

#now construct Gdiff and Omega2
{
    G2_graph <-  graph_from_adjacency_matrix(G2, mode = "undirected")
    G_diff_graph <- igraph::difference(G1_graph, G2_graph)
    G_diff <- as_adjacency_matrix(G_diff_graph)
    {
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
    }

    Omega2 <- Omega1
    #for 50-200 and 100-400
    I_p_sub <- matrix(0, nrow = p_sub, ncol = p_sub)
    diag(I_p_sub) <- diag(Omega2[151:200, 151:200])
    Omega2[151:200, 151:200] <- I_p_sub
    #for 100-200
    #I_p_sub <- matrix(0, nrow = 2 * p_sub, ncol = 2 * p_sub)
    #diag(I_p_sub) <- diag(Omega2[101:200, 101:200])
    #Omega2[101:200, 101:200] <- I_p_sub
    ggcorrplot2(Omega2)
    is.positive.definite(Omega2) 
}

#for 50-400 removing the whole block is too much
#so we'll delete half of edges for each vertex in the last block, starting from highest degree
#trying to keep the subgraph connected and size of the difference to be around 50
{
    G2_subgraph <- G2_subgraph_full <- 
        graph_from_adjacency_matrix(G1[151:200, 151:200],
                                    mode = "undirected")
    gorder(G2_subgraph)
    gsize(G2_subgraph)
    components(G2_subgraph)$no
    
    orig_degrees <- degree(G2_subgraph_full)[order(degree(G2_subgraph_full), 
                                                   decreasing = TRUE)]
    
    for (i in 1:gorder(G2_subgraph)) {
        vertices_by_orig_degree <- V(G2_subgraph)[order(degree(G2_subgraph_full),
                                                        decreasing = TRUE)]
        curr_vertex <- vertices_by_orig_degree[i]
        curr_degree <- degree(G2_subgraph)[curr_vertex]
        if (curr_degree > (orig_degrees[i] / 2)) {
            curr_edges <- incident(G2_subgraph, curr_vertex)
            half_edges <- sample(1:length(curr_edges), 
                                 size = length(curr_edges) - orig_degrees[i] / 2)
            test <- delete_edges(G2_subgraph, 
                                 curr_edges[half_edges])
            test2 <- igraph::difference(G2_subgraph_full, test)
            test2 <- delete_vertices(test2, which(degree(test2) == 0))
            conn_comp_no <- components(test2)$no
            if (conn_comp_no == 1) {
                G2_subgraph <- delete_edges(G2_subgraph, 
                                            curr_edges[half_edges])
                print(gsize(G2_subgraph))
            }
        }
    }
    
    G2_submatrix <- as.matrix(as_adjacency_matrix(G2_subgraph))
    G2 <- G1
    G2[151:200, 151:200] <- G2_submatrix
    
    G2_graph <- graph_from_adjacency_matrix(G2, mode = "undirected")
    G_diff_graph <- igraph::difference(G1_graph, G2_graph)
    gsize(G_diff_graph)
    G_diff <- as.matrix(as_adjacency_matrix(G_diff_graph))
    
    Omega2 <- Omega1 * G2
    diag(Omega2) <- diag(Omega1)
    is.positive.definite(Omega2)
    ggcorrplot2(Omega2)
}


#regardless of the setting, construct Sigma2
Sigma2 <- solve(Omega2)
ggcorrplot2(Sigma2)

save(list = c("G1", "G1_graph", "G2", "G2_graph", 
          "G_diff", "G_diff_graph",
          "Omega1", "Omega2", "Sigma1", "Sigma2",
          "layout"),
     file = "G_Omega_Sigma_scalefree_50_400.Rdata")


