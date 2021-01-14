#!/usr/bin/env Rscript

#########################
### Command line to call the script in linux consol
#########################

"Filter broad peak to keep peaks that passed threshold value

Usage:
  peaks_filter.R [options] <broadPeak_input> <csv_input> <output>
  peaks_filter.R -h | --help

Options:
  -h, --help           Show this screen
" -> doc


library(docopt)
arguments <- docopt(doc)


#########################
### Libraries and input loading
#########################

suppressPackageStartupMessages(suppressWarnings(library(dplyr)))

peaks_bP = read.table(arguments$broadPeak_input, sep = "\t", stringsAsFactors = FALSE)                   # load .broadPeak file
peaks_readcount = read.csv2(arguments$csv_input, sep = ",")                                  # load .csv readcount matrix file


########################
### Script
########################

peaks_to_keep = peaks_readcount  %>% 
  dplyr::filter(nbreads>=10) %>%                                                             # threshold application, suppress peaks that don't pass nbreads >= 10
  dplyr::mutate(peak_threshold_passed = paste(condition,time,donor,peakID,sep="_")) %>%      # new column to have same peak name as in broadPeak file
  dplyr::select(peak_threshold_passed)                                                       # keep only peaks name
peaks_to_keep = peaks_to_keep[[1]]                                                           # transform into character vector

peaks_bP_threshold = peaks_bP %>%                                                            
  dplyr::filter(V4 %in% peaks_to_keep)                                                       # keep lines from original broadPeak that passed the threshold

print("Save .threshold.broadPeak file")  
write.table(peaks_bP_threshold, file = arguments$output, sep ="\t", quote = FALSE)           # save threshold.broadPeak file

