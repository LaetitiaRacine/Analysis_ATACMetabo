
#########################
### Command line to call the script in linux consol
#########################

"Quality control graph

Usage: 
  plot_QC_variation.R -o <png_file> hist_donor <reportcsv_file> <colname>
  plot_QC_variation.R -o <png_file> hist_manip <reportcsv_file> <colname>
  plot_QC_variation.R -o <png_file> line_cond <reportcsv_file> <colname>
  plot_QC_variation.R -o <pdf_file> read_graph <reportcsv_file> <colname> 
  plot_QC_variation.R -o <pdf_file> chrom_single <reportcsv_file> <colname>
  plot_QC_variation.R -o <pdf_file> chrom_multi <reportcsv_file> <colname>
  plot_QC_variation.R -h | --help
  
Options:
  -h, --help               Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)


#########################
### Libraries and input loading
#########################

suppressPackageStartupMessages({
  suppressWarnings(library(ggplot2))
  suppressWarnings(library(viridis))
  suppressWarnings(library(dplyr))
  suppressWarnings(library(ggforce))
  suppressWarnings(library(tidyr))
})

df_report = read.csv2(arguments$reportcsv_file, sep = ";")

#########################
### Histogram plot : facet donor or time, all conditions in one graph
#########################

if (arguments$hist_donor | arguments$hist_manip) {

  var_x = df_report$condition
  var_y = as.numeric(as.character(unlist(df_report[arguments$colname]))) #récupérer la bonne colonne
  var_ylegend = arguments$colname
  
  if (arguments$hist_donor) {
    var_facet = donor~.
    var_fill = df_report$time
    var_legendtitle = "Time spent in culture \n after stimulation"
    var_title = "Inter-donor variability"
  } else if (arguments$hist_manip) {
    var_facet = time~.
    var_fill = df_report$donor
    var_legendtitle = "Donor"
    var_title = "Inter-assay variability "
  }
  
  hist_report = ggplot(df_report, aes(x = var_x, y = var_y, fill = var_fill)) +
    geom_col(color = "black", width = 0.7, position = position_dodge2(width = 0.8, preserve = "single")) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold", colour= "black", size = 16),
          strip.text.x = element_text(size=11, color="black", face="bold.italic"),
          axis.line.y = element_line(colour = "black"),
          axis.text.y = element_text(size = 11, colour = "black"),
          axis.ticks.y = element_line() ,
          axis.title.y = element_text(vjust = 1 ,size = 14),
          axis.title.x = element_blank(),
          axis.text.x = element_text(vjust = -0.5, size = 11, colour = "black"),
          axis.ticks = element_blank()) +
    facet_wrap(var_facet, scales = "free_x") +
    ggtitle(var_title) +
    labs(y = var_ylegend, fill = var_legendtitle) +
    scale_fill_viridis(option="plasma", discrete = TRUE)
  
  ggsave(plot = hist_report,
         filename = arguments$png_file,
         width = 16*0.75, height = 9*0.75)

}

#########################
### Line plot
#########################

if (arguments$line_cond) {

  # On filtre les données pour récupérer seulement celles avec plusieurs points de temps par condition
  df_filtered = data_frame()  
   for (i in 1:length(levels(factor(df_report$condition)))) {
    df_temp = df_report %>%
      dplyr::filter(condition == levels(factor(df_report$condition))[i])
    if (length(levels(factor(df_temp$time)))>1) df_filtered = rbind(df_temp, df_filtered) 
   }
  # On récupère seulement les donneurs avec plusieurs points de temps (facultatif)
  donor_filtered_list = levels(factor(df_filtered$donor))[table(df_filtered$donor)>1]
  df_filtered = df_filtered %>% dplyr::filter(donor %in% donor_filtered_list )
  
  var_x = df_filtered$time
  var_y = as.numeric(as.character(unlist(df_filtered[arguments$colname])))
  var_facet = condition~.
  var_fill = df_filtered$donor
  var_legendtitle = "Donor"
  var_ylegend = arguments$colname
  var_title = "Overtime progression"
  
  line_report = ggplot(df_filtered, aes(x = var_x, y = var_y, group = var_fill)) +
    geom_line(aes(color = var_fill), size = 1)+
    geom_point(aes(color = var_fill), size = 2)+
    theme(legend.position = "right",
          legend.title=element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", colour= "black", size = 16),
          strip.text.x = element_text(size=11, color="black", face="bold.italic"),
          axis.line.y = element_line(colour = "black"),
          axis.text.y = element_text(size = 11, colour = "black"),
          axis.ticks.y = element_line() ,
          axis.title.y = element_text(vjust = 1 ,size = 14),
          axis.title.x = element_blank(),
          axis.text.x = element_text(vjust = -0.5, size = 11, colour = "black"),
          axis.ticks = element_blank()) +
    facet_wrap(var_facet, scales = "free_x") +
    ggtitle(var_title) +
    labs(y = var_ylegend) +
    scale_color_viridis(option="plasma", discrete = TRUE, begin = 0, end = 0.5)
    
  ggsave(plot = line_report,
         filename = arguments$png_file,
         width = 16*0.75, height = 9*0.75)
  
}  


#########################
### Nbreads plot : y peaks count, x nbreads (nbreads_report.csv)
#########################

if (arguments$read_graph) {
  
  df_report = df_report %>%
    dplyr::mutate(condition_time_donor = paste(condition, time, donor, sep="_")) %>%
    dplyr::select(-time, -donor, -peakID) 
  
  var_x = as.numeric(as.character(df_report$nbreads))
  var_ylegend = "Number of peaks"
  
  # Préimpression du graphique pour connaître le nombre de pages à utiliser dans la boucle
  preprint_plot = ggplot(df_report, aes(x=var_x, fill = condition)) +
    geom_histogram(color = "black", position = position_dodge2(width = 0.8, preserve = "single"))+
    geom_vline(xintercept = 10, color ="red", linetype = "dashed", size = 1)+ 
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold", colour= "black", size = 16),
          strip.text.x = element_text(size=11, color="black", face="bold.italic"),
          axis.line.y = element_line(colour = "black"),
          axis.text.y = element_text(size = 11, colour = "black"),
          axis.ticks.y = element_line() ,
          axis.title.y = element_text(vjust = 1 ,size = 14),
          axis.title.x = element_blank(),
          axis.text.x = element_text(vjust = -0.5, size = 11, colour = "black"),
          axis.ticks = element_blank()) +
    xlim(c(0,200)) +
    facet_wrap_paginate(condition_time_donor~., scales = "free_x", ncol = 4, nrow = 4, page = 1) +
    labs(y = var_ylegend) 
  nb_pages = n_pages(preprint_plot)
  
  # Création du vrai graphique
  pdf(width = 16*0.75, height = 9*0.75, file = arguments$pdf_file)
  for(i in 1:nb_pages) {
    print(ggplot(df_report, aes(x=var_x, fill = condition)) +
          geom_histogram(color = "black", position = position_dodge2(width = 0.8, preserve = "single"))+
          geom_vline(xintercept = 10, color ="red", linetype = "dashed", size = 1)+ 
          theme(legend.position = "right",
                plot.title = element_text(hjust = 0.5, face = "bold", colour= "black", size = 16),
                strip.text.x = element_text(size=11, color="black", face="bold.italic"),
                axis.line.y = element_line(colour = "black"),
                axis.text.y = element_text(size = 11, colour = "black"),
                axis.ticks.y = element_line() ,
                axis.title.y = element_text(vjust = 1 ,size = 14),
                axis.title.x = element_blank(),
                axis.text.x = element_text(vjust = -0.5, size = 11, colour = "black"),
                axis.ticks = element_blank()) +
          xlim(c(0,200)) +
          facet_wrap_paginate(condition_time_donor~., scales = "free_x", ncol = 4, nrow = 4, page = i) +
          labs(y = var_ylegend) 
        )
    }
  dev.off()

}


#########################
### Nbpeaks per chromosome plot (peaks_per_chromosome_report.csv)
#########################

if (arguments$chrom_single | arguments$chrom_multi) {

  # Récupération du nom de l'échantillon et du nombre de peak par chromosome
  df_report = df_report %>% 
    dplyr::select(-Start, -End, -Strand) %>% 
    dplyr::mutate(Sample = as.factor(sub("\\_peak.*", "", GeneID))) %>% 
    tidyr::separate(col = Sample, into = c("Condition","Time","Donor"), sep="_", remove = FALSE)
  
  df_recap = data.frame()
  for (i in 1:length(levels(df_report$Sample))) {
    temp = df_report %>% dplyr::filter(Sample == levels(df_report$Sample)[i])
    sample = as.character(temp$Sample[1])
    condition = as.character(temp$Condition[1])
    time = as.character(temp$Time[1])
    donor = as.character(temp$Donor[1])
    
    table = as.data.frame(table(temp$Chr)) %>%
      dplyr::rename("Chr" = Var1, "Peak_count" = Freq) %>%
      dplyr::mutate(Sample = sample, Condition = condition, Time = time, Donor = donor) %>%
      dplyr::mutate(Peak_percentage = round((Peak_count/sum(Peak_count))*100, 2))
    
    df_recap = rbind(df_recap, table)
    rm(temp, table, sample, condition, time, donor)
  }

  # Remplissage de la variable var_y selon le graphique voulu (count/percentage)
  var_y = as.numeric(as.character(unlist(df_recap[arguments$colname])))
  var_ylegend = arguments$colname
  
  
  if (arguments$chrom_single) {   # 1 graphique par chromosome
    
    var_x = df_recap$Sample
    var_fill = df_recap$Sample
    var_facet = Chr~.
          
    # plot  sur 25 pages (car 25 chromosomes)
    pdf(width = 16*0.75, height = 9*0.75, file = arguments$pdf_file)
    for(i in 1:25) {
        print(ggplot(df_recap, aes(x = var_x, y = var_y, fill = var_fill)) +
                geom_col(width = 0.7, position = position_dodge2(width = 0.8, preserve = "single")) +
                theme(legend.position = "none",
                          plot.title = element_text(hjust = 0.5, face = "bold", colour= "black", size = 16),
                          strip.text.x = element_text(size=11, color="black", face="bold.italic"),
                          axis.line.y = element_line(colour = "black"),
                          axis.text.y = element_text(size = 11, colour = "black"),
                          axis.ticks.y = element_line() ,
                          axis.title.y = element_text(vjust = 1 ,size = 14),
                          axis.title.x = element_blank(),
                          axis.text.x = element_text(vjust = -0.5, size = 11, colour = "black", angle = 90),
                          axis.ticks = element_blank()) +
                  labs(y = var_ylegend) +
                  facet_wrap_paginate(var_facet, scales = "free_x", ncol = 1, nrow = 1, page = i) +
                  scale_fill_viridis(option="plasma", discrete = TRUE))
          }
        dev.off()
          
  } else if (arguments$chrom_multi) { # 1 graphique par cond_time avec tous les chromosomes
    
    var_x = df_recap$Chr
    var_fill = df_recap$Donor
    var_facet = Time~Condition
    var_legendtitle = "Donor"
          
    # Preprint pour déterminer le nombre de pages 
    preprint_plot = ggplot(df_recap, aes(x = var_x, y = var_y, fill = var_fill)) +
    geom_col(width = 0.7, position = position_dodge2(width = 0.8, preserve = "single")) +
        theme(legend.position = "right",
              plot.title = element_text(hjust = 0.5, face = "bold", colour= "black", size = 16),
              strip.text.x = element_text(size=11, color="black", face="bold.italic"),
              axis.line.y = element_line(colour = "black"),
              axis.text.y = element_text(size = 11, colour = "black"),
              axis.ticks.y = element_line() ,
              axis.title.y = element_text(vjust = 1 ,size = 14),
              axis.title.x = element_blank(),
              axis.text.x = element_text(vjust = -0.5, size = 9, colour = "black", angle = 90),
              axis.ticks = element_blank()) +
        labs(y = var_ylegend, fill = var_legendtitle) +
        facet_wrap_paginate(var_facet, scales = "free_x", ncol = 2, nrow = 2, page = 1)+
        scale_fill_viridis(option="plasma", discrete = TRUE)
    nb_pages = n_pages(preprint_plot)
          
    # real plot
    pdf(width = 16*0.75, height = 9*0.75, file = arguments$pdf_file)
    for(i in 1:nb_pages) {
        print(ggplot(df_recap, aes(x = var_x, y = var_y, fill = var_fill)) +
                geom_col(width = 0.7, position = position_dodge2(width = 0.8, preserve = "single")) +
                theme(legend.position = "right",
                      plot.title = element_text(hjust = 0.5, face = "bold", colour= "black", size = 16),
                      strip.text.x = element_text(size=11, color="black", face="bold.italic"),
                      axis.line.y = element_line(colour = "black"),
                      axis.text.y = element_text(size = 11, colour = "black"),
                      axis.ticks.y = element_line() ,
                      axis.title.y = element_text(vjust = 1 ,size = 14),
                      axis.title.x = element_blank(),
                      axis.text.x = element_text(vjust = -0.5, size = 9, colour = "black", angle = 90),
                      axis.ticks = element_blank()) +
                labs(y = var_ylegend, fill = var_legendtitle) +
                facet_wrap_paginate(var_facet, scales = "free_x", ncol = 2, nrow = 2, page = i)+
                scale_fill_viridis(option="plasma", discrete = TRUE)
            )
          }
          dev.off()
          
    }
}
