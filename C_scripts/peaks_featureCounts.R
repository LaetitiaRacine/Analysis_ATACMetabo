#!/usr/bin/env Rscript

"Associate reads number to each corresponding region of peaks (number of reads from .bam files)

Usage:
  peaks_featureCounts.R [options] <peaks_gr_input> <bam_input>...
  peaks_featureCounts.R -h | --help

Options:
  --output_rds <file>  Output rds file
  --output_txt <file>  Output matrix count txt file
  --output_csv <file>  Output csv file
  -h, --help           Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)

suppressPackageStartupMessages({
  suppressWarnings(library(GenomicRanges))
  suppressWarnings(library(Rsubread))
  suppressWarnings(library(dplyr))
})

                 
peaks_gr_to_df = function(peaks_gr) {
  peaks_df = data.frame(peaks_gr, stringsAsFactors = FALSE)
  peaks_df$seqnames = as.character(peaks_df$seqnames) # enlever le factor
  peaks_df$strand = as.character(peaks_df$strand)     # enlever le factor
  ranges_df = data.frame(ranges(peaks_gr)) # récuperer le nom des peaks
  peaks_df = full_join(ranges_df, peaks_df, by = c("start","end","width"))
  peaks_df <- peaks_df[,c("names","seqnames","start","end","strand")] # réorganise l'ordre des colonnes
  colnames(peaks_df) = c("GeneID","Chr","Start","End","Strand")
  return(peaks_df)
}

peaks_gr = readRDS(arguments$peaks_gr_input)

peaks_df = peaks_gr_to_df(peaks_gr)


readCount <- featureCounts(files = arguments$bam_input,
                           annot.ext = peaks_df,         # Besoin d'avoir les noms GeneID, Chr, Start, End et Strand pour fonctionner
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

if (!is.null(arguments$output_csv)) {
  print("Save csv...")
  count_df = tibble::rownames_to_column(data.frame(readCount$counts), "name")
  count_df = tidyr::separate(count_df, col = name, into = c("condition","time","donor", "peak", "num"), sep="_", remove = TRUE) 
  count_df = tidyr::unite(count_df, col = peakID, peak, num, sep="_", remove = TRUE)
  colnames(count_df) = c("condition","time","donor","peakID","nbreads")
  write.table(count_df, file = arguments$output_csv, sep = ",", row.names = FALSE) 
  # write.csv2(count_df, file = arguments$output_csv, sep = ",", row.names = FALSE) 
}

