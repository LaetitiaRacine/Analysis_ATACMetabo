#!/usr/bin/env Rscript

"Create and manipulate GRanges files

Usage:
  GRanges.R [options] from_broadPeak <broadPeak_file>
  GRanges.R [options] union <gr_file> <gr_file>...
  GRanges.R [options] intersect <gr_file> <gr_file>...
  GRanges.R -h | --help

Options:
  -o, --output <gr_file>  Output file
  -h, --help              Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)

library(GenomicRanges)
library(dplyr)
library(stringr)

read_peaks_from_macs2 = function(file_path){

  hg19_seqlengths = c("chr1"=249250621, "chr10"=135534747, "chr11"=135006516, "chr12"=133851895, "chr13"=115169878,
                      "chr14"=107349540, "chr15"=102531392, "chr16"=90354753, "chr17"=81195210, "chr18"=78077248,
                      "chr19"=59128983, "chr2"=243199373, "chr20"=63025520, "chr21"=48129895, "chr22"=51304566,
                      "chr3"=198022430, "chr4"=191154276, "chr5"=180915260, "chr6"=171115067, "chr7"=159138663,
                      "chr8"=146364022, "chr9"=141213431, "chrMT"=16571, "chrX"=155270560, "chrY"=59373566)

  # Open .broadPeak file and create columns to further create a Granges file (one per sample)
  peaks <- read.table(file_path,stringsAsFactors = FALSE)
  peaks <- peaks[,c(4,1,2,3)]          # enlève les colonnes inutiles et réorganise l'ordre des colonnes
  colnames(peaks) <- c("GeneID","Chr","Start","End")  # renomme les colonnes
  peaks <- cbind(peaks[,],Strand = rep("*",length(peaks$GeneID))) # ajout de la colonne Strand
  peaks <- data.frame(peaks, stringsAsFactors = FALSE) # évite les problèmes avec les facteurs non voulus

  peaks <- peaks %>%
    mutate(test_chr = str_detect(Chr, pattern = "chr")) %>%
    mutate(test_chrM = str_detect(Chr, pattern = "chrM")) %>%
    mutate(Chr = case_when(test_chrM == TRUE ~ "chrMT", # Transform chrM into chrMT
                           test_chr == TRUE ~ Chr, # Nothing happens
                           test_chr == FALSE ~ paste0("chr", Chr) # If no "chr" are printed in chr column, adds it before chr number
    )
    )

  # Make a GenomicRange object from the data frame "peaks"
  gr_peaks <- GRanges(seqnames = peaks$Chr,
                      ranges = IRanges(peaks$Start, peaks$End),
                      strand = "*")

  names(gr_peaks) = peaks$GeneID # add peak names to Grange (sample_peak_x)
  genome(gr_peaks) <- "hg19" # add genome information to the Grange
  seqlengths(gr_peaks) = hg19_seqlengths # add chr length to each chr

  return(gr_peaks)
}

if (arguments$from_broadPeak) {
  # Transform .broadPeak file into GRanges
  gr = read_peaks_from_macs2(arguments$broadPeak_file)
  # save into output
  saveRDS(gr, file=arguments$output)
}

if (arguments$intersect) {
  gr_list = lapply(X = arguments$gr_file, FUN = readRDS)
  gr_intersection = Reduce(GenomicRanges::intersect, gr_list)
  saveRDS(gr_intersection, file=arguments$output)
}

if (arguments$union) {
  gr_list = lapply(arguments$gr_file, readRDS)
  gr_intersection = Reduce(GenomicRanges::union, gr_list)
  saveRDS(gr_intersection, file=arguments$output)
}
