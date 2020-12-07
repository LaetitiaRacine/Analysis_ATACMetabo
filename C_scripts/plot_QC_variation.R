
#########################
### Création d'une ligne de commande pour appel du script via la console
#########################

"Quality control graph : graph from report nbreads before/after downsampling

Usage:
  plot_QC_variation.R [options] comp_donor <reportcsv_file> <colname>
  plot_QC_variation.R [options] comp_manip <reportcsv_file> <colname>
  plot_QC_variation.R -h | --help
  
Options:
  -o, --output <png_file>  Output file
  -h, --help               Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)


#########################
### Script
#########################

library(ggplot2)

df_report = read.csv2(arguments$reportcsv_file, sep = ",")

var_x = df_report$condition
var_y = unlist(df_report[arguments$colname]) #récupérer la bonne colonne
var_ylegend = arguments$colname

if (arguments$comp_donor) {
  var_facet = donor~.
  var_fill = df_report$time
  var_legendtitle = "Time"
  var_title = "Variabilité interdonneur"
} else if (arguments$comp_manip) {
  var_facet = time~.
  var_fill = df_report$donor
  var_legendtitle = "Donor"
  var_title = "Variabilité intermanip "
}

hist_report = ggplot(df_report, aes(x = var_x, y = var_y, fill = var_fill)) +
  geom_col(color = "black", width = 0.7, position = position_dodge(width = 0.8)) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", colour= "black", size = 16),
        strip.text.x = element_text(size=11, color="black", face="bold.italic"),
        axis.line.y = element_line(colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.ticks.y = element_line() ,
        axis.title.y = element_text(vjust = 1 ,size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = 6, size = 11, colour = "black"),
        axis.ticks = element_blank()) +
  facet_wrap(var_facet, scales = "free_x") +
  ggtitle(var_title) +
  labs(y = var_ylegend, fill = var_legendtitle) +
  scale_fill_viridis_d()

ggsave(plot = hist_report,
       filename = arguments$output,
       width = 10, height = 9)
