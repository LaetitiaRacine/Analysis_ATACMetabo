#########################
### Command line to call the script in linux consol
#########################

"
Usage:
peaks_report.R [options] <input.directory> <output.chrom> <output.glob_long> <output.glob_wide> <output.nbpeaks>
peaks_report.R -h | --help

Options:
-h, --help           Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)

#########################
### Libraries and input loading
#########################

suppressPackageStartupMessages({
  suppressWarnings(library(dplyr))
  suppressWarnings(library(tidyr))
  suppressWarnings(library(stringr))
})


dir.create("D_Analysis/reports")

tab_before = data.frame()
tab_after = data.frame()

dir_readcount = paste0(arguments$input.directory, list.files(arguments$input.directory, pattern = "readcount.csv"))
dir_df_before = paste0(arguments$input.directory, list.files(arguments$input.directory, pattern = "D[[:digit:]].df.csv"))
dir_df_after = paste0(arguments$input.directory, list.files(arguments$input.directory, pattern = "threshold_[[:digit:]]{2,3}.df.csv"))

for (i in 1:length(dir_df_before)) {
  print(paste0("Chargement des fichiers avant threshold : ", i, "/", length(dir_df_before)))
  print(dir_df_before[i])
  df_before = read.csv2(file = dir_df_before[i] , sep = ";") %>%
    tidyr::separate(col = GeneID, into = c("condition", "time", "donor", "peak", "num"), sep = "_", remove = TRUE) %>%
    tidyr::unite(col = sample, condition, time, donor, sep = "_", remove = FALSE) %>%
    tidyr::unite(col = peakID, peak, num, sep = "_", remove = TRUE) %>%
    dplyr::select(-Strand)
  df_before$sample = as.factor(df_before$sample)
  readcount_to_extract = str_extract(string = dir_df_before[i], pattern = "[:alnum:]{2,3}_[:digit:]{2}h_D[:digit:](-D[:digit:])*(?=.df.csv)")
  print(dir_readcount[grep(pattern = readcount_to_extract, dir_readcount)])
  readcount = read.csv2(file = dir_readcount[grep(pattern = readcount_to_extract, dir_readcount)], sep = ";")
  temp = left_join(df_before, readcount, by = c("condition", "time", "donor", "peakID"))
  tab_before = rbind(tab_before,temp)
}
tab_before$threshold = "0"
tab_before$peak = tab_before$peakID

for (i in 1:length(dir_df_after)) {
  print(paste0("Chargement des fichiers avec threshold : ", i, "/", length(dir_df_after)))
  print(dir_df_after[i])
  df_after = read.csv2(file = dir_df_after[i] , sep = ";") %>%
    tidyr::separate(col = GeneID, into = c("condition", "time", "donor", "peak", "num"), sep = "_", remove = TRUE) %>%
    tidyr::unite(col = sample, condition, time, donor, sep = "_", remove = FALSE) %>%
    tidyr::unite(col = peakID, peak, num, sep = "_", remove = TRUE) %>%
    dplyr::select(-Strand) %>%
    dplyr::mutate(threshold = str_extract(dir_df_after[i], pattern = "(?<=threshold_)[:digit:]{2,3}(?=.df.csv)"))
  df_after$sample = as.factor(df_after$sample)
  readcount_to_extract = str_extract(string = dir_df_after[i], pattern = "[:alnum:]{2,3}_[:digit:]{2}h_D[:digit:](-D[:digit:])*(?=_threshold_[:digit:]{2,3}.df.csv)")
  print(dir_readcount[grep(pattern = readcount_to_extract, dir_readcount)])
  readcount = read.csv2(file = dir_readcount[grep(pattern = readcount_to_extract, dir_readcount)], sep = ";")
  temp = full_join(df_after,
                    readcount %>% filter(peakID %in% df_after$peakID),
                    by = c("condition", "time", "donor", "peakID"))
  tab_after = rbind(tab_after, temp)
}

rm(df_before, df_after, readcount, temp)

########################
### Data frame creation
########################

## Tableaux propres pour lecture

tab_global_wide = full_join(tab_before %>% tidyr::pivot_wider(names_from = threshold, names_prefix = "threshold_", values_from = peakID), 
                            tab_after %>% tidyr::pivot_wider(names_from = threshold, names_prefix = "threshold_", values_from = peakID),
                            by = c("sample","condition", "time", "donor", "Chr", "Start", "End", "nbreads")) %>% 
                  mutate(across(.cols = contains("threshold"), .fns = ~!is.na(.)))


tab_global_long =  full_join(tab_before, 
                             tab_after,
                             by = c("sample","condition", "time", "donor", "peakID", "nbreads", "Chr", "Start", "End", "threshold")) %>%
                   dplyr::select(-peak)
tab_global_long$threshold = as.factor(tab_global_long$threshold)
  

## Tableau nbpeaks_report pour graphiques

tab_nbpeaks = tab_global_long %>%
  group_by(sample, condition, time, donor, threshold) %>%
  summarise(nbpeaks = n(), mean_nbreads = mean(nbreads)) 

initial_value = tab_nbpeaks %>% slice_min(threshold)
tab_nbpeaks$initial = rep(initial_value$nbpeaks, each = length(levels(tab_nbpeaks$threshold)))
tab_nbpeaks = tab_nbpeaks %>% 
  mutate(lost_percentage = ifelse(threshold == "0", 0, round(((initial-nbpeaks)/initial)*100, 4))) %>%
  select(-initial)

rm(initial_value)

## Tableau number_peak_per_chromosome_report

tab_chromosome = tab_global_long %>%
  group_by(sample, condition, time, donor, threshold, Chr) %>%
  summarise(nbpeaks = n())

########################
### Output saving
########################

write.table(tab_chromosome , file = arguments$output.chrom, sep = ";", row.names = FALSE) 
write.table(tab_global_long , file = arguments$output.glob_long, sep = ";", row.names = FALSE) 
write.table(tab_global_wide , file = arguments$output.glob_wide, sep = ";", row.names = FALSE) 
write.table(tab_nbpeaks , file = arguments$output.nbpeaks, sep = ";", row.names = FALSE) 
