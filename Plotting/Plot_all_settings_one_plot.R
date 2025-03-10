library(ggplot2)

perf$both_sizes <- paste(perf$true_diff_size, perf$true_G1_size,
                         sep = "-")
perf$both_sizes <- factor(perf$both_sizes,
                          levels = c("50-200", "100-200",
                                     "50-400", "100-400"))

perf$sample_size <- factor(perf$sample_size, 
                           levels = c("400 samples",
                                      "100 samples"))


text_size <- 26
title_size <- 24
linewidth <- 1.2

ggplot(perf[perf$method == "DTrace", ], 
       aes(x = FDR, y = power, 
           linetype = both_sizes, color = graph_type)) +
    geom_line(linewidth = linewidth) + 
    facet_grid(~sample_size) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = c(0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1)) +
    ylim(0, 1) +
    scale_linetype_manual(values = c(1, 2, 3, 5)) +
    theme(text = element_text(size = text_size),
          plot.title = element_text(size = title_size),
          legend.key.size = unit(3.5, "lines")) 
