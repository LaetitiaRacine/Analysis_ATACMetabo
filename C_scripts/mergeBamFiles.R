#########################
### Command line to call the script in linux consol
#########################

" Merge bam from different donors but same condition and time

Usage: 
  mergeBamFiles.R --output_dir <merged_dir> <csv_file> <bam_dir>
 
Options:
  -h, --help               Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)

#########################
### Libraries and input loading
#########################

suppressPackageStartupMessages({
  suppressWarnings(library(dplyr))
  suppressWarnings(library(stringr))
  suppressWarnings(library(Rsamtools))
})

df_linktab = read.csv2(arguments$csv_file, sep = ",")

#########################
### Script
#########################

df_linktab = df_linktab %>% 
  dplyr::filter(USE == "true") %>%
  mutate(condition_time = paste0(CONDITION, "_", TIME))

cond_time = levels(as.factor(df_linktab$condition_time))
  
for (i in 1:length(cond_time)) {
  # Keep line corresponding to one condition_time
      df_temp = df_linktab %>% dplyr::filter(condition_time == cond_time[i])
  # Extract donors name
      donors = paste(levels(factor(df_temp$DONOR)),collapse = "-")
  # Extract input files names to merge
      file_names = paste0(arguments$bam_dir, "/", levels(factor(df_temp$RAW.FILE)))
      print(file_names)
  # Create output name 
      output_name = paste0(arguments$merged_dir, cond_time[i], "_", donors, ".bam")
      print(output_name)
  # Merge input bam files and save the result in chosen output_name
      mergeBam(file_names, output_name)
  }


