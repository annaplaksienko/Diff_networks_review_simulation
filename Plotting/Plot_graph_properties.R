library(igraph)
library(ggplot2)
library(latex2exp)

load("Generate_data/G1_G2_Gdiff_Sigma1_Sigma2_all.Rdata")
load("Plotting/perf_AUC.Rdata")

sapply(G_diff_list, diameter)
sapply(G_diff_list, radius)
sapply(G_diff_list, girth)
sapply(G_diff_list, global_efficiency)
sapply(G_diff_list, average_local_efficiency)
sapply(G_diff_list, count_components)
sapply(G_diff_list, function(x) max(degree(x)))
sapply(G_diff_list, clique_num)
sapply(G_diff_list, count_motifs)
sapply(G_diff_list, function(x) authority_score(x)$value)
sapply(G_diff_list, function(x) hub_score(x)$value)
sapply(G_diff_list, mean_distance)

p <- 200
perf_AUC$value <- sapply(G_diff_list, 
                          function(x) max(degree(x)) / p)

text_size <- 22

graph_type_Tex_title <- TeX("$|G^{diff}|-|G^{(1)}|$")  

#perf_AUC <- perf_AUC %>% arrange(desc(graph_type))

ggplot(perf_AUC, aes(x = value, y = AUC, 
                     color = both_sizes, shape = graph_type)) +
    geom_point(size = 6, alpha = 0.9) +
    ylim(c(0, 1)) +
    theme(text = element_text(size = text_size)) +
    xlab("Highest degree / p") +
    ylab("Mean AUC") +
    scale_shape_manual(name = "Graph type",
                       values = c(15, 17, 8)) +
    scale_color_discrete(name = graph_type_Tex_title)
