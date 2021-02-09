#### Supplementary figure 1 - ANGSD Het GATK (-GL 2) vs Samtools (-GL 1) ####

# Load packages
library(ggplot2)
library(reshape2)

# Load data into r
He_indBuf <- read.csv("19-Heterozygosity_buffalo.csv", header = T)
# Reformat data to better work with ggplot2
He_reshape <- melt(He_indBuf, value.name = "Heterozygosity")

# Plot heterozygosity comparison per sample:
Het_GL2vGL1 <- ggplot(He_reshape) +
  geom_line(aes(x=Sample, y=Heterozygosity), col = "grey") +
  geom_point(aes(x=Sample, y=Heterozygosity, fill = variable), pch = 21, size = 2) +
  scale_x_discrete(limits=He_indBuf$Sample) + #Plot in the same order as in data file, and not alphabetically
  scale_fill_manual(name = "Genotype likelihood model:",
                    labels = c("-GL 1 (Samtools)", "-GL 2 (GATK)"),
                    values = c("#d7191c", "#2c7bb6")) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        legend.position = "top",
        legend.spacing.x = unit(0.05, "cm")) +
  labs(y="Genome-wide heterozygosity") +
  coord_flip() #Since the input file is in reverse alphabetical order, Addo samples will be at the top after coord_flip
Het_GL2vGL1

# Save as pdf
pdf("Het_GL2vGL1.pdf",
    width = 11.69,
    height = 8.27)
Het_GL2vGL1
dev.off()
