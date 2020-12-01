#!/usr/bin/env Rscript

"Add annotations in GRanges files

Usage:
  annotate_grange.R [--output_csv <csv_file>] -a <annotations_file> <grange_input_file> <grange_annotated_file>
  annotate_grange.R -h | --help

Options:
  --output_csv <csv_file>  Output annotated file
  -h, --help               Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)

# # Test in RStudio (without using console)
# ATAC_dir="/path/to/Analysis_ATACMetabo/"
# gr_input_dir=paste0(ATAC_dir, "D_results/genomic_ranges/static_peaks/")
# output_dir=paste0(ATAC_dir, "D_results/genomic_ranges/static_peaks_annotated/")
# arguments <- docopt(doc, args=c(
#   "--output_csv", paste0(output_dir, "AOA_24h_D1,D2,D3_ann.gr.csv"),
#   "-a", paste0(ATAC_dir, "A_raw_data/Annotation_TSS_pm1kb_int_ex_53utr_ctcf_cpg_histo_gr.rda"),
#   paste0(gr_input_dir, "AOA_24h_D1,D2,D3.gr.rds"),
#   paste0(output_dir, "AOA_24h_D1,D2,D3_ann.gr.rds")
# ))

library(GenomicRanges)
library(dplyr)
# library(stringr)

loadRData<-function(fileName){
  #loads an RData file, and returns it  
  load(fileName)
  get(ls()[ls()!="fileName"])
}

# the "big" file containing "all" the annotations 
all_annotations = loadRData(arguments$annotations_file)

gr = readRDS(arguments$grange_input_file)

# First a matrix is created filled with FALSE and added to the Grange
annotations_types = levels(factor(all_annotations$annotation))
metadata = matrix(FALSE, ncol = length(annotations_types), nrow = length(gr))
colnames(metadata) = annotations_types
mcols(gr) = metadata # ajout de metadata : une colonne par annotation toutes à false au démarrage

# for each of the annotations types an overlap is calculated and used to assigned the peak as TRUE when overlapping with the annotation
for (i in 1:ncol(metadata)){
  sub_annot = all_annotations[all_annotations$annotation == annotations_types[i]]
  overlaps = findOverlaps(gr, sub_annot)
  mcols(gr)[queryHits(overlaps),i] = TRUE
}

colnames(mcols(gr)) = c("UTR3P","UTR5P","CpG", "CTCF","Exons","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3","Introns","TSS_mp1kb")

mcols(gr) = as_tibble(mcols(gr)) %>%
  dplyr::mutate(Intergenic = ifelse(UTR3P == FALSE & UTR5P == FALSE & Exons == FALSE & Introns == FALSE & TSS_mp1kb == FALSE, TRUE, FALSE)) %>%
  dplyr::mutate(CpG_Intergenic = ifelse(Intergenic == TRUE & CpG == TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(CpG_Intergenic = ifelse(Intergenic == TRUE & CpG == TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(CTCF_Intergenic = ifelse(Intergenic == TRUE & CTCF == TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(CTCF_in_intron = ifelse(Introns == TRUE & CTCF == TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(CTCF_in_exon = ifelse(Exons == TRUE & CTCF == TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(Histone_Intergenic = ifelse(Intergenic== TRUE & (H3K27me3 | H3K36me3 | H3K4me1 | H3K4me3 | H3K9me3) == TRUE, TRUE, FALSE))

saveRDS(gr, arguments$grange_annotated_file)
if (!is.null(arguments$output_csv)) {
  write.csv2(as.data.frame(gr), arguments$output_csv)
}
