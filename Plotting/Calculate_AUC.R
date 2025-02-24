library(ggplot2)
library(dplyr)

#remember to set the path to the folder with the code and the data
#path <- getwd()
load(paste(path, "/Plotting/perf_ALL.Rdata", sep = ""))

perf$both_sizes <- paste(perf$true_diff_size, perf$true_G1_size,
                         sep = "-")
perf$both_sizes <- factor(perf$both_sizes,
                          levels = c("50-200", "100-200",
                                     "50-400", "100-400"))
perf <- perf[perf$sample_size == "400 samples", ]

perf_AUC <- perf %>% 
    group_by(graph_type, both_sizes, sample_size) %>%
    summarise(AUC = MESS::auc(FDR, power))
perf_AUC

save(list = c('perf_AUC'), file = "perf_AUC.Rdata")
