load("Generate_data/G1_G2_Gdiff_Sigma1_Sigma2_all.Rdata")
load("Plotting/perf_AUC.Rdata")

distance <- rep(NA, 12)
names(distance) <- names(G_diff_list)

library(dad)
library(latex2exp)

for (i in 1:12) {
    #library(ProDenICA)
    #distance[i] <- (amari(Sigma1_list[[i]], Sigma2_list[[i]]) + amari(Sigma2_list[[i]], Sigma1_list[[i]])) / 2
    
    #library(fpc)
    p <- 200
    #distance[i] <- bhattacharyya.dist(rep(0, p), rep(0, p), Sigma1_list[[i]], Sigma2_list[[i]])
    
    #library(dad)
    #distance[i] <- wassersteinpar(rep(0, p), Sigma1_list[[i]], rep(0, p), Sigma2_list[[i]])
    #distance[i] <- hellingerpar(rep(0, p), Sigma1_list[[i]], rep(0, p), Sigma2_list[[i]])
    distance[i] <- jeffreyspar(rep(0, p), Sigma1_list[[i]], rep(0, p), Sigma2_list[[i]])
    #distance[i] <- l2dpar(rep(0, p), Sigma1_list[[i]], rep(0, p), Sigma2_list[[i]], check = TRUE)
}

perf_AUC$distance <- distance

text_size <- 22

graph_type_Tex_title <- TeX("$|G^{diff}|-|G^{(1)}|$")  

ggplot(perf_AUC, aes(x = log2(distance), y = AUC, 
                 color = both_sizes, shape = graph_type)) +
    geom_point(size = 6, alpha = 0.9) +
    ylim(c(0, 1)) +
    theme(text = element_text(size = text_size)) +
    xlab("Log2 of symm. KL divergence") +
    ylab("Mean AUC") +
    scale_shape_manual(name = "Graph type",
                       values = c(15, 17, 8)) +
    scale_color_discrete(name = graph_type_Tex_title)
