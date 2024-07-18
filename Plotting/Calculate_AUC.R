library(ggplot2)
library(dplyr)

load("perf_ALL.Rdata")

perf <- rbind(perf[, 2:13], perf_FGL[4:15])
perf$both_sizes <- paste(perf$true_diff_size, perf$true_G1_size,
                         sep = "-")
perf$both_sizes <- factor(perf$both_sizes,
                          levels = c("50-200", "100-200",
                                     "50-400", "100-400"))
perf <- perf[perf$sample_size == "400 samples", ]
perf <- perf[perf$TP + perf$FP > 1, ]

ggplot(perf, aes(x = FDR, y = power, 
                 color = both_sizes,
                 linetype = graph_type)) +
    geom_line(linewidth = 1.5) +
    facet_grid(~method) +
    xlim(0, 1) + ylim(0, 1) 

ggplot(perf, aes(x = FDR, y = power, 
                 color = method)) +
    geom_line(linewidth = 1.5) +
    facet_grid(graph_type ~ both_sizes) +
    xlim(0, 1) + ylim(0, 1) 

perf_AUC <- perf %>% 
    group_by(graph_type, both_sizes, sample_size) %>%
    summarise(AUC = MESS::auc(FDR, power))

save(list = c('perf_AUC'), file = "perf_AUC.Rdata")
