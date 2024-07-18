library(igraph)
#for is.positive.definite
library(matrixcalc)

source("ggcorrplot_with_limits.R")

#number of vertices
p <- 200
#number of subnetworks
L <- 4
p_sub <- p / L

m <- 2
power <- 1

vsize <- 2
par(mar = c(1, 1, 1, 1) + 0.5)

Omega_sub <- vector(mode = "list", length = L)
for (l in 1:L) {
    G <- as.matrix(as_adjacency_matrix(barabasi.game(p_sub, 
                                                     power = power, m = m, 
                                                     directed = FALSE)))
    size <- sum(G) / 2
    
    a <- 0.6
    b <-  0.9
    samp_right <- runif(ceiling(size / 2), min = a, max = b)
    samp_left <- runif(ceiling(size / 2), min = -b, max = -a)
    samp <- sample(c(samp_left, samp_right), size)
    ind <- 1
    for (i in 2:p_sub) {
        for (j in 1:(i - 1)) {
            if (G[i, j] != 0) {
                G[i, j] <- G[j, i] <- samp[ind]
                ind <- ind + 1
            }
        }
    }
    
    Omega_sub[[l]] <- G
}

Omega1 <- matrix(0, nrow = p, ncol = p)
for (l in 1:L) {
    Omega1[((l - 1) * p_sub + 1):(l * p_sub), 
           ((l - 1) * p_sub + 1):(l * p_sub)] <- Omega_sub[[l]]
}
ggcorrplot2(Omega1, limits = range(Omega1))

G1 <- Omega1 != 0
G1_graph <- graph_from_adjacency_matrix(G1, mode = "undirected")
gsize(G1_graph)
layout <- layout_nicely(G1_graph)
plot(G1_graph, vertex.frame.color = "white", 
     layout = layout,
     vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
     edge.width = 3)

is.positive.definite(Omega1)
diag(Omega1) <- rowSums(abs(Omega1)) + 0.1
is.positive.definite(Omega1)
ggcorrplot2(Omega1, limits = range(Omega1))

Sigma1 <- solve(Omega1)
ggcorrplot2(Sigma1, limits = range(Omega1))

save(list = c("G1", "G1_graph", "Omega1", "Sigma1", "layout"),
     file = "G1_Omega1_Sigma1_m2_power1_scalefree_random.Rdata")
