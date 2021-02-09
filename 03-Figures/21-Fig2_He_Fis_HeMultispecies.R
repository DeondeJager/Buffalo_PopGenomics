#### Fig 2 - Multiple plots on one page with ggplot2 and cowplot ####

# Install the following packages
#Install cowplot
#install.packages("cowplot")
#install.packages("svglite")

# Load packages
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggrepel)
library(svglite)

# Load data into r
He_indBuf <- read.csv("19-Heterozygosity_buffalo.csv", header = T)
Fis <- read.csv("21-Inbreeding_coefficients.csv", header = T)
He_species <- read.csv("21-Heterozygosity_other_species&buffalo.csv", header=T)

## Make plot objects

# Specify the factor levels for IUCNstatus in the order you want (otherwise will be alphabetical)
He_species$IUCNstatus <- factor(He_species$IUCNstatus, levels = c("EW", "CR", "EN", "VU", "NT", "LC", "DD"))

# Scatterplot of species heterozygosity (Part C of Fig2) 
# Draw the plot and save it as an object ("Ryan used me as an object" - Kelly Kapoor)
# Used https://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=5 to get 5 middle colours for the palette
He_species_point <- ggplot(He_species, aes(x=Census, y=Heterozygosity)) +
  geom_pointrange(aes(ymin=Min_Het, ymax=Max_Het), size = 0.5) +
  geom_text_repel(aes(x = Census, y = Heterozygosity, label = Common.name), point.padding = 0.1) +
  geom_point(aes(x = Census, y = Heterozygosity, fill = IUCNstatus), size = 3, pch = 21) +
  theme_classic() +
  scale_fill_manual(values = c('black','#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6', 'white'),
                    name = "IUCN") +
  scale_x_log10(breaks = c(1e1, 1e2, 1e3, 1e4, 1e5), labels = c("10", "100", "1,000", "10,000", "100,000")) +
  theme(text = element_text(size = 12)) +
  labs(x="Census size", y="Genome-wide heterozygosity")
He_species_point


# Heterozygosity of individual buffalo (Part A of Fig2)
He_buf_plot <- ggplot(He_indBuf, aes(Sample, Het_GL1, fill = Pop)) +
  scale_fill_manual(values = c("#0066FF","#CC0000","#666666","#FF3300"),
                    breaks = c("AENP", "KNP", "HiP", "MNP")) +
  geom_col(colour = "black", size = 0.3) +
  scale_x_discrete(limits=He_indBuf$Sample) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.0035)) +
  theme_classic() +
  theme(text = element_text(size=12),
        axis.text.x = element_text(),
        axis.title.x = element_text(),
        axis.title.y = element_text(vjust = 2),
        axis.text.y = element_text(size = 6), 
        legend.position = "None") +
  labs(y="Genome-wide heterozygosity") +
  coord_flip()
He_buf_plot

# Inbreeding coefficient of individual buffalo (Part B of Fig2)
Fis_buf_plot <- ggplot(Fis, aes(Sample, Inbreeding.coefficient, fill = Pop)) +
  scale_fill_manual(values = c("#0066FF","#CC0000","#666666","#FF3300"),
                    breaks = c("AENP", "KNP", "HiP", "MNP"),
                    name = "Population") +
  geom_col(colour = "black", size = 0.3) +
  scale_x_discrete(limits=Fis$Sample) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(text = element_text(size=12),
        axis.text.x = element_text(),
        axis.title.x = element_text(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  labs(y=expression(paste("Individual inbreeding coefficient (", italic("F"),")"))) +
  coord_flip()
Fis_buf_plot


# Now use the plot_grid() function of cowplot to arrange plots
# See: https://wilkelab.org/cowplot/articles/plot_grid.html

# Align plots in grid
plots <- align_plots(He_buf_plot, Fis_buf_plot, align = "h", axis = "t", greedy = FALSE)

top_panel <- plot_grid(plots[[1]], plots[[2]],
                        labels = c("a", "b"), hjust = c(-0.5, 0.8),
                        nrow = 1, ncol = 2)
top_panel
final_plot <- plot_grid(top_panel, He_species_point,
                        labels = c("", "c"), 
                        nrow = 2)
final_plot

# Save as svg file
# Then opened in Inkscape to bold Cape buffalo label and move text labels around for clarity.
ggsave2("He_Fis_HeMultispecies_new_ggsave2.svg", #ggsave2 does not use Dingbats font
       plot = final_plot,
       width = 11.69,
       height = 8.27,
       units = "in")

# Can also save directly as a pdf, without bold text for Cape buffalo and some labels may obscure some points
#pdf("He_Fis_HeMultispecies_new.pdf",
#    width = 11.69,
#    height = 8.27)
#final_plot
#dev.off()
