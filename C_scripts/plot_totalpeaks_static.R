#!/usr/bin/env Rscript

"Create plot exploiting static peaks annotated grange (total nb peaks per condition)

Usage:
  plot_totalpeaks_static.R [options] <output_plot> <annotated_gr>...
  plot_totalpeaks_static.R -h | --help

Options:
  -h, --help           Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)

#################################################################
################### Libraries ###################################
#################################################################

library(GenomicRanges)
library(ggplot2)
library(ggthemes)
library(stringr)
library(dplyr)

list_gr = lapply(arguments$annotated_gr, readRDS)

time_condition_from_filename = function(file_path) {
  name = basename(file_path)
  name_split = unlist(str_split(name, "_"))
  condition = name_split[1]
  time = name_split[2]
  return(paste(time, condition, sep=" "))
}

df_nb_peaks = tibble("nb_peaks" = unlist(lapply(list_gr, length)),
                     "time_condition" = unlist(lapply(arguments$annotated_gr, time_condition_from_filename)))

total_peak_hist <- ggplot(df_nb_peaks, aes(x = time_condition, y = nb_peaks)) +
  geom_col(aes(fill = time_condition), color = "black") +
  theme_tufte()+
  theme(legend.position = "none",
        axis.line.y = element_line(colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.ticks.y = element_line() ,
        axis.title.y = element_text(vjust = 1 ,size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = 6, size = 11, colour = "black"),
        axis.ticks = element_blank())+
  #ggtitle(paste0("Culture condition ", conditions_list[drug])) +
  geom_text(aes(label = nb_peaks), position = position_dodge(width = 0.8), vjust = -0.25)+
  ylim(c(0,max(df_nb_peaks$nb_peaks) + 5000)) +
  labs(y = "Total number of peaks detected") +
  scale_fill_viridis_d()
total_peak_hist

ggsave(plot = total_peak_hist, filename = arguments$output_plot, width = 10, height = 9)
