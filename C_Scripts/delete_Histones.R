# Script pour enlever les histones du Grange d'annotations 

#***********************
# Chargement des fonctions et librairies
#***********************

library(dplyr)

loadRData<-function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls()!="fileName"])
}

#**********************
# Chargement du Grange initial
#**********************

initial = as.data.frame(loadRData(fileName = "Annotation_TSS_pm1kb_int_ex_53utr_ctcf_cpg_histo_gr.rda"))

#**********************
# Effacement des lignes correspondants aux histones
#**********************

woh = initial %>% dplyr::filter(annotation != "H3K9me3",
                                annotation != "H3K36me3",
                                annotation != "H3K27me3", 
                                annotation != "H3K4me3",
                                annotation != "H3K4me1")


hg19_seqlengths = c("chr1"=249250621, "chr2"=243199373,"chr3"=198022430, "chr4"=191154276, "chr5"=180915260, "chr6"=171115067,
                    "chr7"=159138663, "chr8"=146364022, "chr9"=141213431, "chr10"=135534747, "chr11"=135006516, "chr12"=133851895, 
                    "chr13"=115169878, "chr14"=107349540, "chr15"=102531392, "chr16"=90354753, "chr17"=81195210, "chr18"=78077248,
                    "chr19"=59128983, "chr20"=63025520, "chr21"=48129895, "chr22"=51304566, "chrMT"=16571, "chrX"=155270560, "chrY"=59373566)
chromosomes_represented = unique(woh$seqnames)
hg19_seqlengths = hg19_seqlengths[names(hg19_seqlengths) %in% chromosomes_represented]

gr = GRanges(seqnames = woh$seqnames,
             ranges = IRanges(woh$start, woh$end),
             strand = woh$strand,
             annotation = woh$annotation, 
             genome = "hg19",
             seqlengths = hg19_seqlengths)

saveRDS(gr, file="Annotation_TSS_pm1kb_int_ex_53utr_ctcf_cpg_woThisto_gr.rda")
