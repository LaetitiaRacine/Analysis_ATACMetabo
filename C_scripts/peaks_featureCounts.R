#!/usr/bin/env Rscript

"Associate reads number to each corresponding region of peaks (number of reads from .bam files)

Usage:
  peaks_featureCounts.R [options] <peaks_gr_input> <bam_input>...
  peaks_featureCounts.R -h | --help

Options:
  --output_rds <file>  Output rds file
  --output_txt <file>  Output matrix count txt file
  -h, --help           Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)

library(Rsubread)

peaks_gr_to_df = function(gr) {
  peaks_df = data.frame(peaks_gr, stringsAsFactors = FALSE)
  peaks_df$seqnames = as.character(peaks_df$seqnames) # enlever le factor
  peaks_df$strand = as.character(peaks_df$strand)     # enlever le factor
  peaks_df$name = paste0("peak_",rownames(peaks_df))  # ajout colonne nom du peak
  peaks_df <- peaks_df[,c("name","seqnames","start","end","strand")] # réorganise l'ordre des colonnes
  colnames(peaks_df) <- c("GeneID","Chr","Start","End","Strand") # renomme les colonnes
  return(peaks_df)
}

peaks_gr = readRDS(arguments$peaks_gr_input)

peaks_df = peaks_gr_to_df(peaks_gr)

readCount <- featureCounts(files = arguments$bam_input,
                           annot.ext = peaks_df,
                           isPairedEnd = TRUE,
                           nthreads = 1,
                           countChimericFragments = FALSE,
                           countMultiMappingReads = TRUE
)

if (!is.null(arguments$output_rds)) {
  print("Save rds...")
  saveRDS(readCount, file = arguments$output_rds)
}

if (!is.null(arguments$output_txt)) {
  print("Save txt...")
  matrix_count = readCount$counts
  write.table(matrix_count, file = arguments$output_txt, sep = "\t", quote = FALSE)
}