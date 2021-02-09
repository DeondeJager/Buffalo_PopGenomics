#### Fig 4 - Main PSMC plot ####
# Script adapted for this data set from Emily Humble's ("3.1_psmc.R" from https://github.com/elhumble/SHO_analysis_2020)
# Calls the "psmc.result" function from the "plot_psmc.R" script from https://figshare.com/articles/Plot_PSMC_results/3996156/1 

#~~ Load requried packages
library(ggplot2)
source("plot_psmc.R") # Load plot_psmc function. File must be in current working directory, else give path to file.
library(data.table)
library(plyr)
library(tidyr)
library(scales)
options(scipen=999) # Disable scientific notation of large numbers- see "Scipen" under ?options

#~~ Specify variables for plotting
i.iteration=25 # Number of iterations to use in file
s=100 # bin size

#~~ Read in main PSMC files. 
psmc_files <- paste("data/psmc/", list.files(path = "data/psmc", pattern="*.psmc"), sep = "")

#~~ Set mu and g (for plot_psmc function):
mu <- 1.5e-8
g <- 7.5

#~~ Run "psmc.result" from the "plot_psmc" function
psmc_buf <- lapply(psmc_files, psmc.result, i.iteration, mu = mu, s, g = g)

#~~ Set names in psmc_buf to sample names
dataset_names <-list.files(path = "data/psmc/", pattern="*.psmc")
dataset_names <-  gsub(".psmc", "", dataset_names)
names(psmc_buf) <- dataset_names

#~~ Transform list into dataframe
psmc_buf <- ldply(psmc_buf, .id = "Sample")  

#~~ Get bootstraps
boot_files <- paste("data/psmc_boot/", list.files(path = "data/psmc_boot", pattern="*round*"), sep = "")

#~~ Run "psmc.result" from the "plot_psmc" function
boot <- lapply(boot_files, psmc.result, i.iteration, mu, s, g)

#~~ Set names in psmc_buf to sample names
dataset_names <-list.files(path = "data/psmc_boot/", pattern="*round*")
dataset_names <-  gsub(".psmc", "", dataset_names)
names(boot) <- dataset_names

#~~ Transform list into dataframe
boot_df <- ldply(boot, .id = "ID") %>%
  separate(ID, c("Sample", "Boot"), sep = "_round-", remove = F) %>%
  mutate(Boot = gsub("_", "", Boot))

#~~ Define colours
cbPalette <- c("#0066FF", "#CC0000", "#666666", "#FF3300")

#~~ Make plot and store as object
psmc_main <- ggplot() +
  geom_line(aes(YearsAgo, Ne, group = Boot), 
            dplyr::filter(boot_df, Sample == "A_243_14" & mu == mu & g == g), 
            size = .1, alpha = 0.2, col = cbPalette[1]) +
  geom_line(aes(YearsAgo, Ne, color = Sample),  
            dplyr::filter(psmc_buf, Sample == "A_243_14"), size = 0.7, alpha=0.9) +
  geom_line(aes(YearsAgo, Ne, group = Boot), 
            dplyr::filter(boot_df, Sample == "B98_509" & mu == mu & g == g), 
            size = .1, alpha = 0.2, col = cbPalette[2]) +
  geom_line(aes(YearsAgo, Ne, color = Sample), 
            dplyr::filter(psmc_buf, Sample == "B98_509"), size = 0.7, alpha=0.9) +
  geom_line(aes(YearsAgo, Ne, group = Boot), 
            dplyr::filter(boot_df, Sample == "HC_32" & mu == mu & g == g), 
            size = .1, alpha = 0.2, col = cbPalette[3]) +
  geom_line(aes(YearsAgo, Ne, color = Sample),
            dplyr::filter(psmc_buf, Sample == "HC_32"), size = 0.7, alpha=0.9) +
  geom_line(aes(YearsAgo, Ne, group = Boot), 
            dplyr::filter(boot_df, Sample == "M_120_13" & mu == mu & g == g), 
            size = .1, alpha = 0.2, col = cbPalette[4]) +
  geom_line(aes(YearsAgo, Ne, color = Sample),
            dplyr::filter(psmc_buf, Sample == "M_120_13"), size = 0.7, alpha=0.9) +
  theme_classic() +
  scale_colour_manual(values = cbPalette,
                      labels = c("AENP", "KNP", "HiP", "MNP"),
                      name = expression(paste(bold("Population")))) +
  theme(legend.position = c(0.95,0.8)) +
  scale_x_log10(label = comma,
                breaks = c(10000,20000,50000,100000,200000,500000,1000000,2000000)) +
  annotation_logticks(sides = "b", alpha = 0.5) +
  scale_y_continuous(label = comma) +
  xlab(expression(paste("Years before present (g = 7.5, ", mu," = 1.5e-8",")"))) + 
  ylab (expression(paste("IICR (scaled in units of 4",italic(N[e]),mu,")"))) +
  geom_vline(xintercept = 193000, linetype = "dashed", colour = "grey") + # Subspecies split
  annotate("text", x = 175000, y=275000, label = "Subspecies split", angle = 90, size = 5) +
  #annotate("rect", xmin = 11700, xmax = 115000, ymin = 310000, ymax = Inf,
           #alpha = .1) + # Last glacial period
  annotate("rect", xmin = 75000, xmax = 135000, ymin = -Inf, ymax = Inf, 
           alpha = 0.1) + # Droughts
  annotate("text", x=100000, y=275000, label = "Mega-droughts", angle = 90, size =5) +
  geom_vline(xintercept=22000, linetype = "dashed", colour = "grey") + # LGM
  #annotate("text", x = 20000, y = 275000, label = "LGM 22 ka", angle = 90, size = 5) +
  annotate("rect", xmin=12000, xmax=48000, ymin = -Inf, ymax = Inf,
           alpha = 0.1) + # Colonisation of southern Africa
  annotate("text", x=23000, y=310000, label="Colonisation of\nsouthern Africa", size =5) +
  geom_vline(xintercept=11700, colour = "grey") + # Holocene
  annotate("text", x = 10900, y = 275000, label = "Holocene", angle = 90, size = 5) +
  geom_vline(xintercept=2580000, colour = "grey") + # Pleistocene
  annotate("text", x = 2380000, y = 275000, label = "Pleistocene", angle = 90, size = 5)
psmc_main

#~~ Save as pdf:
pdf(file="figs/psmc_figure.pdf", 
    height=8.27, width=11.69)
psmc_main
dev.off()

#~~ Save as svg:
ggsave("figs/psmc_figure.svg",
        plot = psmc_main,
        width = 11.69,
        height = 8.27,
        units = "in")
  