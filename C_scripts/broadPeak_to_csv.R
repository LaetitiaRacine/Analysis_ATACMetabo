
#**********************
# Command line to call the script in linux consol
#**********************

"Transform .broadPeak file into organized data.frame (.csv) to be used with FeatureCount or ready to be transformed int oGrange

Usage:
broadPeak_to_csv.R [options] <broadPeak_input> <output>
broadPeak_to_csv.R -h | --help

Options:
-h, --help           Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)

#**********************
# Libraries and data loading
#**********************

suppressPackageStartupMessages({
  suppressWarnings(library(dplyr))
  suppressWarnings(library(stringr))
})

broadPeak_file = read.table(arguments$broadPeak_input, sep ="\t", stringsAsFactors = FALSE)

#**********************
# Script to create df
#**********************

broad_df <- broadPeak_file[,c(4,1,2,3)] 
colnames(broad_df) <- c("GeneID","Chr","Start","End")
broad_df <- cbind(broad_df[,],Strand = rep("*",length(broad_df$GeneID)))
broad_df <- data.frame(broad_df)   
broad_df$GeneID = as.character(broad_df$GeneID) # enlève les factors
broad_df$Chr = as.character(broad_df$Chr) # enlève les factors
broad_df <- broad_df %>%
  dplyr::mutate(test_chr = str_detect(Chr, pattern = "chr")) %>%
  dplyr::mutate(test_chrM = str_detect(Chr, pattern = "chrM")) %>%
  dplyr::mutate(Chr = case_when(test_chrM == TRUE ~ "chrMT", # Transform chrM into chrMT
                         test_chr == TRUE ~ Chr, # Nothing happens
                         test_chr == FALSE ~ paste0("chr", Chr)))  %>% # If no "chr" are printed in chr column, adds it before chr number
  dplyr::select(-test_chr, -test_chrM)      

write.table(broad_df, file = arguments$output, sep = ",", row.names = FALSE) 
