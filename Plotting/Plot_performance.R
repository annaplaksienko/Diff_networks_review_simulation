library(ggplot2)
library(latex2exp)

load("perf_ALL.Rdata")

perf_random <- perf[perf$graph_type == "random", ]
perf_random_100 <- perf_random[perf_random$sample_size == "100 samples", ]
perf_random_400 <- perf_random[perf_random$sample_size == "400 samples", ]
perf_sf <- perf[perf$graph_type == "scalefree", ]
perf_sf_100 <- perf_sf[perf_sf$sample_size == "100 samples", ]
perf_sf_400 <- perf_sf[perf_sf$sample_size == "400 samples", ]
perf_star <- perf[perf$graph_type == "star", ]
perf_star_100 <- perf_star[perf_star$sample_size == "100 samples", ]
perf_star_400 <- perf_star[perf_star$sample_size == "400 samples", ]

perf_random_FGL <- perf_FGL[perf_FGL$graph_type == "random", ]
perf_random_100_FGL <- perf_random_FGL[perf_random_FGL$sample_size == "100 samples", ]
perf_random_400_FGL <- perf_random_FGL[perf_random_FGL$sample_size == "400 samples", ]
perf_sf_FGL <- perf_FGL[perf_FGL$graph_type == "scalefree", ]
perf_sf_100_FGL <- perf_sf_FGL[perf_sf_FGL$sample_size == "100 samples", ]
perf_sf_400_FGL <- perf_sf_FGL[perf_sf_FGL$sample_size == "400 samples", ]
perf_star_FGL <- perf_FGL[perf_FGL$graph_type == "star", ]
perf_star_100_FGL <- perf_star_FGL[perf_star_FGL$sample_size == "100 samples", ]
perf_star_400_FGL <- perf_star_FGL[perf_star_FGL$sample_size == "400 samples", ]

appender1 <- function(string) {
    TeX(paste("$|G^{diff}| \\approx $", string), output = "character")  
}

appender2 <- function(string) {
    TeX(paste("$|G^{(1)}| \\approx $", string), output = "character")  
}

diff_labels_unique <- unique(appender1(perf$true_diff_size))
diff_labels <- c("50" = diff_labels_unique[2], 
                 "100" = diff_labels_unique[1])

text_size <- 28
title_size <- 24
linewidth <- 1.2

ggplot(perf_random_100, aes(x = FDR, y = power, color = method)) +
    geom_line(linewidth = linewidth) + 
    facet_grid(appender2(true_G1_size) ~ 
                   factor(appender1(true_diff_size),
                          levels = diff_labels), 
               labeller = labeller(.rows = label_parsed, 
                                   .cols = label_parsed)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = c(0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1)) +
    ylim(0, 1) +
    scale_color_manual(values = c("#00BA38", "#00BFC4", 
                                  "#619CFF", "#F564E3", "gold", 
                                  "#F8766D")) +
    theme(text = element_text(size = text_size),
          plot.title = element_text(size = title_size),
          legend.position = "none") +
    geom_path(data = perf_random_100_FGL, 
              aes(x = FDR, y = power), linewidth = linewidth)

ggplot(perf_random_400, aes(x = FDR, y = power, color = method)) +
    geom_line(linewidth = linewidth) + 
    facet_grid(appender2(true_G1_size) ~ 
                   factor(appender1(true_diff_size),
                          levels = diff_labels), 
               labeller = labeller(.rows = label_parsed, 
                                   .cols = label_parsed)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = c(0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1)) +
    ylim(0, 1) +
    scale_color_manual(values = c("#00BA38", "#00BFC4", 
                                  "#619CFF", "#F564E3", "gold", 
                                  "#F8766D")) +
    theme(text = element_text(size = text_size),
          plot.title = element_text(size = title_size),
          legend.position = "none") +
    geom_path(data = perf_random_400_FGL, 
              aes(x = FDR, y = power), linewidth = linewidth)

ggplot(perf_sf_100, aes(x = FDR, y = power, color = method)) +
    geom_line(linewidth = linewidth) + 
    facet_grid(appender2(true_G1_size) ~ 
                   factor(appender1(true_diff_size),
                          levels = diff_labels), 
               labeller = labeller(.rows = label_parsed, 
                                   .cols = label_parsed)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = c(0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1)) +
    ylim(0, 1) +
    scale_color_manual(values = c("#00BA38", "#00BFC4", 
                                  "#619CFF", "#F564E3", "gold", 
                                  "#F8766D")) +
    theme(text = element_text(size = text_size),
          plot.title = element_text(size = title_size),
          legend.position = "none") +
    geom_path(data = perf_sf_100_FGL, 
              aes(x = FDR, y = power), linewidth = linewidth)

ggplot(perf_sf_400, aes(x = FDR, y = power, color = method)) +
    geom_line(linewidth = linewidth) + 
    facet_grid(appender2(true_G1_size) ~ 
                   factor(appender1(true_diff_size),
                          levels = diff_labels), 
               labeller = labeller(.rows = label_parsed, 
                                   .cols = label_parsed)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = c(0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1)) +
    ylim(0, 1) +
    scale_color_manual(values = c("#00BA38", "#00BFC4", 
                                  "#619CFF", "#F564E3", "gold", 
                                  "#F8766D")) +
    theme(text = element_text(size = text_size),
          plot.title = element_text(size = title_size),
          legend.position = "none") +
    geom_path(data = perf_sf_400_FGL, 
              aes(x = FDR, y = power), linewidth = linewidth)

ggplot(perf_star_100, aes(x = FDR, y = power, color = method)) +
    geom_line(linewidth = linewidth) + 
    facet_grid(appender2(true_G1_size) ~ 
                   factor(appender1(true_diff_size),
                          levels = diff_labels), 
               labeller = labeller(.rows = label_parsed, 
                                   .cols = label_parsed)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = c(0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1)) +
    ylim(0, 1) +
    scale_color_manual(values = c("#00BA38", "#00BFC4", 
                                  "#619CFF", "#F564E3", "gold", 
                                  "#F8766D")) +
    theme(text = element_text(size = text_size),
          plot.title = element_text(size = title_size),
          legend.position = "none") +
    geom_path(data = perf_star_100_FGL, 
              aes(x = FDR, y = power), linewidth = linewidth)

ggplot(perf_star_400, aes(x = FDR, y = power, color = method)) +
    geom_line(linewidth = linewidth) + 
    facet_grid(appender2(true_G1_size) ~ 
                   factor(appender1(true_diff_size),
                          levels = diff_labels), 
               labeller = labeller(.rows = label_parsed, 
                                   .cols = label_parsed)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = c(0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1)) +
    ylim(0, 1) +
    scale_color_manual(values = c("#00BA38", "#00BFC4", 
                                  "#619CFF", "#F564E3", "gold", 
                                  "#F8766D")) +
    theme(text = element_text(size = text_size),
          plot.title = element_text(size = title_size),
          legend.position = "none") +
    geom_path(data = perf_star_400_FGL, 
              aes(x = FDR, y = power), linewidth = linewidth) 
