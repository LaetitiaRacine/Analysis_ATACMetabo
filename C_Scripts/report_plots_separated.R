
#########################
### Command line to call the script in linux consol
#########################

"Reports value in graph

Usage: 

  report_plots_separated.R -o <output_name> <input_report_file> <graph_type> <extension> <colname>
  report_plots_separated.R -h | --help
  
Options:
  -h, --help               Show this screen

" -> doc

library(docopt)
arguments <- docopt(doc)


#########################
### Libraries, functions and input loading
#########################

suppressPackageStartupMessages({
  suppressWarnings(library(ggplot2))
  suppressWarnings(library(viridis))
  suppressWarnings(library(dplyr))
  suppressWarnings(library(ggforce))
  suppressWarnings(library(tidyr))
  suppressWarnings(library(stringr))
})

df_report = read.csv2(arguments$input_report_file, sep = ";")

############################################################
### Graphs from qc_report et nbreads_report
############################################################

if (arguments$input_report_file == "D_Analysis/reports/qc_report.csv" | arguments$input_report_file == "D_Analysis/reports/nbreads_report.csv") {

  ####### Code spécifique pour ajouter la valeur de downsampling en sous-titre
  
  if (arguments$colname == "nbreads_after_downsampling") { 
    var_subtitle = paste("Downsampling value = ", df_report$downsampling_value[1]) 
  } else { 
    var_subtitle = "" 
  }
  
  ####### Code spécifique pour le tracé des histogrammes 
  
  if (arguments$graph_type == "hist_time" | arguments$graph_type == "hist_donor") {
    
    var_x = df_report$condition
    var_y = as.numeric(as.character(unlist(df_report[arguments$colname])))
    var_ylegend = arguments$colname

    if (arguments$graph_type == "hist_time") {
      var_facet = time~.
      var_fill = df_report$donor
      var_legendtitle = "Donor"
      var_title = "Inter-timing variability" 
    } else {
      var_facet = donor~.
      var_fill = df_report$time
      var_legendtitle = "Time spent in culture \n after stimulation"
      var_title = "Inter-donor variability"
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
      labs(title = var_title, subtitle = var_subtitle, y = var_ylegend, fill = var_legendtitle)+
      scale_fill_viridis(option="plasma", discrete = TRUE)
    
    ggsave(plot = hist_report,
           filename = arguments$output_name,
           width = 16*0.75, height = 9*0.75)
    
  }
  
  ####### Code spécifique pour le tracé des lignes de points
  
  if (arguments$graph_type == "line_cond" | arguments$graph_type == "line_donor") {
    
    # On filtre les données pour récupérer seulement celles avec plusieurs points de temps par condition
    df_filtered = data_frame()  
    for (i in 1:length(unique(df_report$condition))) {
      df_temp = df_report %>% dplyr::filter(condition == unique(df_report$condition)[i])
      if (length(unique(df_temp$time))>1) df_filtered = rbind(df_temp, df_filtered) 
    }
    rm(df_temp)
    
    # On attribue les valeurs dans les variables
    var_x = df_filtered$time
    var_y = as.numeric(as.character(unlist(df_filtered[arguments$colname])))
    var_title = "Overtime progression"
    var_ylegend = arguments$colname
    
    if (arguments$graph_type == "line_cond") {
      var_facet = condition~.
      var_fill = df_filtered$donor
      var_legendtitle = "Donor"
    } else { 
      var_facet = donor~.
      var_fill = df_filtered$condition
      var_legendtitle = "Condition"
    }
    
    line_report = ggplot(df_filtered, aes(x = var_x, y = var_y, group = var_fill)) +
      geom_line(aes(color = var_fill), size = 1)+
      geom_point(aes(color = var_fill), size = 2)+
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
      labs(title = var_title, subtitle = var_subtitle, y = var_ylegend)+
      scale_color_viridis(option="plasma", discrete = TRUE, begin = 0, end = 0.5, name = var_legendtitle)
    
    ggsave(plot = line_report,
           filename = arguments$output_name,
           width = 16*0.75, height = 9*0.75)
  }
  
}

############################################################
### Graphs from nbpeaks_per_chromosome_report 
############################################################

if (arguments$input_report_file == "D_Analysis/reports/nbpeaks_per_chromosome_report.csv") {
  
  if (arguments$extension == "png") {
    
    print("Impossible de créer un graphique png sur plusieurs pages : changer l'extension en .pdf.")
    
  } else {
    
    for (i in 1:length(unique(df_report$threshold))) {

      print(paste("threshold", unique(df_report$threshold)[i]))
      
      df_threshold = df_report  %>% 
        dplyr::rename("peak_count" = nbpeaks) %>%
        dplyr::filter(threshold == unique(df_report$threshold)[i])  %>%
        group_by(sample, condition, time, donor, threshold) %>%
        dplyr::mutate(sum_peaks = sum(peak_count)) %>%
        dplyr::mutate(peak_percentage = round((peak_count/sum_peaks)*100, 2)) %>%
        dplyr::select(-sum_peaks)
      
      # Extraction of saving file's name
      name = str_extract(arguments$output_name, pattern = paste0(".+(?=_threshold_[0-9]+.", arguments$extension, ")"))
      extension = str_sub(arguments$output_name, nchar(arguments$output_name)-3, nchar(arguments$output_name))
      threshold = paste0("threshold_", unique(df_threshold$threshold))
      
      var_subtitle = paste0("Threshold = ", unique(df_threshold$threshold))
      var_y = as.numeric(as.character(unlist(df_threshold[arguments$colname])))
      var_ylegend = arguments$colname
      
      if (arguments$graph_type == "chrom_single") {   # 1 graphique par chromosome
        
        var_x = df_threshold$sample
        var_fill = df_threshold$sample
        var_facet = Chr~.
        
        # plot  sur 25 pages (car 25 chromosomes)
        pdf(width = 16*0.75, height = 9*0.75, file = paste0(name, "_", threshold, extension))
        for(j in 1:25) {   # changer le 25 en 1 quand seulement chromosomeX
          print(ggplot(df_threshold, aes(x = var_x, y = var_y, fill = var_fill)) +
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
                  labs(y = var_ylegend, subtitle = var_subtitle) +
                  facet_wrap_paginate(var_facet, scales = "free_x", ncol = 1, nrow = 1, page = j) +
                  scale_fill_viridis(option="plasma", discrete = TRUE))
        }
        dev.off()
        
      } else if (arguments$graph_type == "chrom_multi") { # 1 graphique par cond_time avec tous les chromosomes
        
        var_x = df_threshold$Chr
        var_fill = df_threshold$donor
        var_facet = time~condition
        var_legendtitle = "Donor"
        
        # Preprint pour déterminer le nombre de pages 
        preprint_plot = ggplot(df_threshold, aes(x = var_x, y = var_y, fill = var_fill)) +
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
          labs(y = var_ylegend, fill = var_legendtitle, subtitle = var_subtitle) +
          facet_wrap_paginate(var_facet, scales = "free_x", ncol = 2, nrow = 2, page = 1)+
          scale_fill_viridis(option="plasma", discrete = TRUE)
        nb_pages = n_pages(preprint_plot)
        
        # real plot
        pdf(width = 16*0.75, height = 9*0.75, file = paste0(name, "_", threshold, extension))
        for(j in 1:nb_pages) {
          print(ggplot(df_threshold, aes(x = var_x, y = var_y, fill = var_fill)) +
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
                  labs(y = var_ylegend, fill = var_legendtitle, subtitle = var_subtitle) +
                  facet_wrap_paginate(var_facet, scales = "free_x", ncol = 2, nrow = 2, page = j)+
                  scale_fill_viridis(option="plasma", discrete = TRUE)
          )
        }
        dev.off()
      }
    }
  }
}

############################################################
### Graphs from nbpeaks_report 
############################################################

if (arguments$input_report_file == "D_Analysis/reports/nbpeaks_report.csv") {
  
  if (arguments$extension == "png") {
    
    print("Impossible de créer un graphique png sur plusieurs pages : changer l'extension en .pdf.")
    
  } else {
    
    df_report$mean_nbreads = as.numeric(df_report$mean_nbreads)
    df_report$lost_percentage = as.numeric(df_report$lost_percentage)
    var_subtitle = ""
    
    if (length(unique(df_report$threshold))>4) { 
      npage = round(length(unique(df_report$threshold))/4) 
    } else { 
      npage = 1 
    }
    
    if (arguments$graph_type == "hist_time" | arguments$graph_type == "hist_donor") {
      
      var_x = df_report$condition
      var_y = as.numeric(as.character(unlist(df_report[arguments$colname])))
      var_ylegend = arguments$colname

      if (arguments$graph_type == "hist_time") {
        col = length(unique(df_report$time))
        var_facet = threshold~time
        var_fill = df_report$donor
        var_legendtitle = "Donor"
        var_title = "Inter-timing variability"
      } else {
        col = length(unique(df_report$donor))
        var_facet = threshold~donor
        var_fill = df_report$time
        var_legendtitle = "Time spent in culture \n after stimulation"
        var_title = "Inter-donor variability"
      }

      pdf(width = 16*0.75, height = 9*0.75, file = arguments$output_name)
      for(i in 1:npage) {
        print(ggplot(df_report, aes(x = var_x, y = var_y, fill = var_fill)) +
                geom_col(color = "black", width = 0.7, position = position_dodge2(width = 0.8, preserve = "single")) +
                theme(legend.position = "right",
                      plot.title = element_text(hjust = 0.5, face = "bold", colour= "black", size = 16),
                      strip.text.x = element_text(size=11, color="black", face="bold.italic"),
                      axis.line.y = element_line(colour = "black"),
                      axis.text.y = element_text(size = 11, colour = "black"),
                      axis.ticks.y = element_line() ,
                      axis.title.y = element_text(vjust = 1 ,size = 14),
                      axis.title.x = element_blank(),
                      axis.text.x = element_text(vjust = -0.5, size = 9, colour = "black"),
                      axis.ticks = element_blank()) +
                labs(title = var_title, subtitle = var_subtitle, y = var_ylegend, fill = var_legendtitle)+
                scale_fill_viridis(option="plasma", discrete = TRUE) +
                facet_grid_paginate(var_facet, scales = "free_x", page = i, nrow = 4, ncol = col))
                # pas possible de mettre des variables pour nrow et ncol sinon erreur dans snakemake
      }
      dev.off()
    }
    
    if (arguments$graph_type == "line_cond" | arguments$graph_type == "line_donor") {

      # On filtre les données pour récupérer seulement celles avec plusieurs points de temps par condition
      df_filtered = data_frame()  
      for (i in 1:length(unique(df_report$condition))) {
        df_temp = df_report %>% dplyr::filter(condition == unique(df_report$condition)[i])
        if (length(unique(df_temp$time))>1) df_filtered = rbind(df_temp, df_filtered) 
      }
      rm(df_temp)

      # On attribue les valeurs dans les variables
      var_x = df_filtered$time
      var_y = as.numeric(as.character(unlist(df_filtered[arguments$colname])))
      var_title = "Overtime progression"
      var_ylegend = arguments$colname

      if (arguments$graph_type == "line_cond") {
        col = length(unique(df_filtered$condition))
        var_facet = threshold~condition
        var_fill = df_filtered$donor
        var_legendtitle = "Donor"
      } else {
        col = length(unique(df_filtered$donor))
        var_facet = threshold~donor
        var_fill = df_filtered$condition
        var_legendtitle = "Condition"
      }
  
      pdf(width = 16*0.75, height = 9*0.75, file = arguments$output_name)
      for(i in 1:npage) {
        print(ggplot(df_filtered, aes(x = var_x, y = var_y, group = var_fill)) +
                geom_line(aes(color = var_fill), size = 1)+
                geom_point(aes(color = var_fill), size = 2)+
                theme(legend.position = "right",
                      plot.title = element_text(hjust = 0.5, face = "bold", colour= "black", size = 16),
                      strip.text.x = element_text(size=11, color="black", face="bold.italic"),
                      axis.line.y = element_line(colour = "black"),
                      axis.text.y = element_text(size = 11, colour = "black"),
                      axis.ticks.y = element_line() ,
                      axis.title.y = element_text(vjust = 1 ,size = 14),
                      axis.title.x = element_blank(),
                      axis.text.x = element_text(vjust = -0.5, size = 11, colour = "black", angle = 90),
                      axis.ticks = element_blank()) +
                facet_grid_paginate(var_facet, scales = "free_x", page = i, nrow = 4, ncol = col) +
                labs(title = var_title, subtitle = var_subtitle, y = var_ylegend)+
                scale_color_viridis(option="plasma", discrete = TRUE, begin = 0, end = 0.5, name = var_legendtitle))
      }
      dev.off()
    }
  }
}

#########################
### Graphs from nbpeaks_nbreads_report
#########################

if (arguments$input_report_file == "D_Analysis/reports/nbpeaks_nbreads_long_report.csv") {
  
  if (arguments$extension == "png") {
    
    print("Impossible de créer un graphique png sur plusieurs pages : changer l'extension en .pdf.")
    
  } else {
    
    var_x = as.numeric(as.character(df_report$nbreads))
    var_ylegend = "Number of peaks"
    
    # Préimpression du graphique pour connaître le nombre de pages à utiliser dans la boucle
    preprint_plot = ggplot(df_report, aes(x=var_x, fill = condition)) +
      geom_histogram(color = "black", position = position_dodge2(width = 0.8, preserve = "single"))+
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
      facet_wrap_paginate(sample~., scales = "free_x", ncol = 4, nrow = 4, page = 1) +
      labs(y = var_ylegend) 
    nb_pages = n_pages(preprint_plot)
    
    # Création du vrai graphique
    pdf(width = 16*0.75, height = 9*0.75, file = arguments$output_name)
    for(i in 1:nb_pages) {
      print(ggplot(df_report, aes(x=var_x, fill = condition)) +
              geom_histogram(color = "black", position = position_dodge2(width = 0.8, preserve = "single"))+
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
              facet_wrap_paginate(sample~., scales = "free_x", ncol = 4, nrow = 4, page = i) +
              labs(y = var_ylegend))
    }
    dev.off()
  }
}    

  
 








