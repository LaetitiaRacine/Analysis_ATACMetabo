#!/usr/bin/env Rscript

"Differential analysis between two datasets
(DEseq2, foldChange and pValue associated calculation)

Usage:
  DEseq.R [options] <readCount_file>
  DEseq.R -h | --help

Options:
  -d, --deseqresult <file>  Output DESeq result file
  -h, --help                Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)

readCount = readRDS(arguments$readCount_file)
# Only keep number of reads detected for each samples in the union region
matrix_count = readCount$counts


#library(GenomicRanges)
library(DESeq2)
#library(Rsubread)
library(dplyr)
#library(stringr)

### Parameters  for DEseq
seqType <- rep("paired-end", ncol(matrix_count))
condition <- c(rep("condition1", (ncol(matrix_count)/2)), rep("condition2", (ncol(matrix_count)/2)))
# on part du principe que pour comparer deux time_point, on doit avoir le même nombre de sample dans chaque
coldata <- as.data.frame(cbind(condition = condition,  type = seqType), row.names = colnames(matrix_count))
# Prepare DEseq object
dds <- DESeqDataSetFromMatrix(matrix_count, colData = coldata, design = ~ condition)
# Filter peaks base on minimum reads
keep <- rowSums(counts(dds) >= 10) >= 2 # At least 2 samples must have a minimum of 10 reads to consider the region
dds <- dds[keep,]
# Calculate DEseq Matrix
dds <- DESeq(dds)
res <- results(dds)



# Merging peaks informations and DEseq results
readCount_annot = readCount$annotation %>%
  dplyr::rename(name = GeneID) %>%
  dplyr::select(-Length,  -Strand)
deseq_annot <- as.data.frame(res)
deseq_annot$name <- rownames(deseq_annot)
deseq_annot <- merge(readCount_annot, deseq_annot, by=c("name")) # Adding chr/start/end informations
if (sum(is.na(deseq_annot$pvalue)) != 0) deseq_annot = deseq_annot[-which(is.na(deseq_annot$pvalue)),]   ### Remove rows where DEseq generated NA for pValues
## Indicate differential status of the region between 2 conditions
deseq_annot = deseq_annot %>%
  mutate(regulation = case_when(log2FoldChange < 0 ~ "gained-close",
                                log2FoldChange > 0 ~ "gained-open"))

# Export final results
if (!is.null(arguments$deseqresult)) {
  write.table(deseq_annot, file = arguments$deseqresult, sep = "\t", quote = FALSE)
}


############
### 4: Create genomic.ranges from DEseq2 results
############

Peak_Grange = makeGRangesFromDataFrame(deseq_annot %>% select (name, Chr, Start, End, regulation),  # Keeping only peak info + regulation
                                       keep.extra.columns = TRUE,
                                       ignore.strand = FALSE,
                                       seqinfo = NULL,
                                       seqnames.field= c ("seqnames", "seqname",
                                                          "chromosome", "chrom",
                                                          "chr", "chromosome_name",
                                                          "seqid"),
                                       start.field = "start",
                                       end.field = c("end", "stop"),
                                       strand.field = "strand",
                                       starts.in.df.are.0based = FALSE)

# Create and save increasing regions
peaks_inc_gr = Peak_Grange[Peak_Grange$regulation == "gained-open"]
peaks_inc_gr = clean_gr(peaks_inc_gr)

save(peaks_inc_gr,
     file = paste0(dir_genomic_ranges_differential,"peaks_inc_",
                   name_condition_list_b, "_vs_", name_condition_list_a, "_gr.rda"))

# Create and save decreasing regions
peaks_dec_gr = Peak_Grange[Peak_Grange$regulation == "gained-close"]
peaks_dec_gr = clean_gr(peaks_dec_gr)


save(peaks_dec_gr,
     file = paste0(dir_genomic_ranges_differential,"peaks_dec_",
                   name_condition_list_b, "_vs_", name_condition_list_a, "_gr.rda"))
