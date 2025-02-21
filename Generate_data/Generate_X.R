#for mvrnorm
library(MASS)

load("G1_G2_Gdiff_Sigma1_Sigma2_all.Rdata")

L <- length(G_diff_list)

X_list <- vector(mode = "list", length = L)
names(X_list) <- names(G_diff_list)

p <- 200
nreps <- 50
n_vec <- c(100, 400)

for (l in 1:L) {
    print(l)
    X_list[[l]] <- vector(mode = "list", length = 2)
    names(X_list[[l]]) <- paste("n", n_vec, sep = "")
    
    Sigma1 <- Sigma1_list[[l]]
    Sigma2 <- Sigma2_list[[l]]
    
    for(n in c(1:length(n_vec))) {
        X_reps <- vector(mode = "list", length = nreps)
        names(X_reps) <- paste("Replication", 1:nreps, sep = "")
        for (r in 1:nreps) {
            X_reps[[r]] <- vector(mode = "list", length = 2)
            names(X_reps[[r]]) <- c("X1", "X2")
            X_reps[[r]][[1]] <- mvrnorm(n = n_vec[n], mu = rep(0, p), Sigma = Sigma1)
            X_reps[[r]][[2]] <- mvrnorm(n = n_vec[n], mu = rep(0, p), Sigma = Sigma2)
        }
        X_list[[l]][[n]] <- X_reps
    }
    
    remove(Sigma1, Sigma2, X_reps)
}

save(list = c("G_diff_list", 'X_list', "settings"),
     file = "G_diff_X_all.Rdata")
