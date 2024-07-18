library(igraph)
#for is.positive.definite
library(matrixcalc)

source("ggcorrplot_with_limits.R")
#load("G1_Omega1_Sigma1_m2_400_power1_scalefree_random.Rdata")

p <-  200
#L <-  4
#p_sub <-  p / L

G2_graph <- rewire(G1_graph, keeping_degseq(niter = 30))
G_diff_graph <- igraph::union(igraph::difference(G1_graph, G2_graph), 
                        igraph::difference(G2_graph, G1_graph))
gsize(G_diff_graph)
G2 <- as.matrix(as_adjacency_matrix(G2_graph))
G_diff <- as.matrix(as_adjacency_matrix(G_diff_graph))
{
    common <- intersection(G1_graph, G2_graph)
    E(common)$color <- 1
    
    diff21 <- add_edges(common, 
                        as.vector(t(as_edgelist(difference(G2_graph, G1_graph)))),
                        color = "2")
    diff12 <- add_edges(common, 
                        as.vector(t(as_edgelist(difference(G1_graph, G2_graph)))),
                        color = "2")
    vsize <- 2
    
    plot(diff12, vertex.frame.color = "white", 
         layout = layout, 
         vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
         edge.width = 3)
    
    plot(diff21, vertex.frame.color = "white", 
         layout = layout, 
         vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
         edge.width = 3)
    
    E(G_diff_graph)$color <- 2
    plot(G_diff_graph, vertex.frame.color = "white", 
         layout = layout, 
         vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
         edge.width = 3)
    
}

ggcorrplot2(Omega1, limits = range(Omega1))
diag(G2) <- 1
Omega2 <- Omega1 * G2
ggcorrplot2(Omega2, limits = range(Omega2))
diag(G2) <- 0
ggcorrplot2(G2, limits = range(G2))

size <- gsize(G_diff_graph)

a <- 0.6
b <- 0.9
samp_right <- runif(ceiling(size / 2), min = a, max = b)
samp_left <- runif(ceiling(size / 2), min = -b, max = -a)
samp <- sample(c(samp_left, samp_right), size)
ind <-  1
for (i in 2:p) {
    for (j in 1:(i - 1)) {
        if (G2[i, j] != 0 & G1[i, j] == 0) {
            Omega2[i, j] <- Omega2[j, i] <- samp[ind]
            ind <- ind + 1
        }
    }
}

is.positive.definite(Omega1)
is.positive.definite(Omega2)
#diag(Omega2) <- diag(Omega2) + 0.1
is.positive.definite(Omega2)
ggcorrplot2(Omega1, limits = range(Omega1))
ggcorrplot2(Omega2, limits = range(Omega2))

ggcorrplot2(Sigma1, limits = range(Sigma1))
Sigma2 <-  solve(Omega2)
ggcorrplot2(Sigma2, limits = range(Sigma2))


save(list = c("G1", "G1_graph", "G2", "G2_graph", 
          "G_diff", "G_diff_graph",
          "Omega1", "Omega2", "Sigma1", "Sigma2",
          "layout"),
     file = "G_Omega_Sigma_random_100_400.Rdata")


