#### Supplementary figure 2 - Heterozygosity vs mean coverage ####

# Load packages
library(ggplot2)
library(reshape2)
library(dplyr)

# Load data into r
Het_cov <- read.csv("20-Subsampled_heterozygosity.csv", header = T)

# Reorder coverage category for legend
Het_cov$CoverageCategory <- factor(Het_cov$CoverageCategory, 
                                   levels = c("Low (0.25 reads)", "Medium (0.5 reads)", "High (1.0 reads)"))

# Plot
Het_vs_meanCov <- ggplot(Het_cov) +
  geom_point(aes(x=Mean.coverage, y=Heterozygosity, fill = CoverageCategory), pch = 21, size = 3) +
  facet_wrap(~Sample) +
  labs(x="Mean coverage", y="Genome-wide heterozygosity") +
  scale_fill_manual(values = c("#d7191c", "#ffffbf", "#2c7bb6"),
                    name = "Coverage category\n(proportion reads subsampled):",
                    labels = c("Low (25% reads)", "Medium (50% reads)", "High (no subsampling)")) +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "top",
        legend.title.align = 0.5)

# Save as pdf
pdf("Het_vs_meanCov.pdf",
    width = 11.69,
    height = 8.27)
Het_vs_meanCov
dev.off()
