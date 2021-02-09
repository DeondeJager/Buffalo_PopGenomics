#### Supplementary figure 5 - Individual PSMC plots, coloured by population ####
# Script adapted for this data set from Emily Humble's ("3.1_psmc.R" from https://github.com/elhumble/SHO_analysis_2020)
# Calls the "psmc.result" function from the "plot_psmc.R" script from https://figshare.com/articles/Plot_PSMC_results/3996156/1 

#~~ Load required packages
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
psmc_files <- paste("data/psmc_all/", list.files(path = "data/psmc_all", pattern="*.psmc"), sep = "")

#~~ Set mu and g (for plot_psmc function):
mu <- 1.5e-8
g <- 7.5

#~~ Run "psmc.result" from the "plot_psmc" function
psmc_buf <- lapply(psmc_files, psmc.result, i.iteration, mu = mu, s, g = g)

#~~ Set names in psmc_buf to sample names
dataset_names <-list.files(path = "data/psmc_all/", pattern="*.psmc")
dataset_names <-  gsub(".psmc", "", dataset_names)
names(psmc_buf) <- dataset_names

#~~ Transform list into dataframe
psmc_buf <- ldply(psmc_buf, .id = "Sample")  

#~~ Add population info by checking the first letter of the sample name and assigning it to the given pop if true and
# store as a new dataframe "psmc_df"
psmc_df <- psmc_buf %>%
  mutate(pop = ifelse(grepl("^A", Sample), "AENP",
                       ifelse(grepl("^B", Sample), "KNP",
                              ifelse(grepl("^H", Sample), "HiP",
                                     ifelse(grepl("^M", Sample), "MNP", "NA")))))


#~~ Get bootstraps
boot_files <- paste("data/psmc_boot_all/", list.files(path = "data/psmc_boot_all", pattern="*round*"), sep = "")

#~~ Run "psmc.result" from the "plot_psmc" function
boot <- lapply(boot_files, psmc.result, i.iteration, mu, s, g)

#~~ Set names in psmc_buf to sample names
dataset_names <-list.files(path = "data/psmc_boot_all/", pattern="*round*")
dataset_names <-  gsub(".psmc", "", dataset_names)
names(boot) <- dataset_names

#~~ Transform list into dataframe and add pop info
boot_df <- ldply(boot, .id = "ID") %>%
  separate(ID, c("Sample", "Boot"), sep = "_round-", remove = F) %>%
  mutate(Boot = gsub("_", "", Boot)) %>%
  mutate(pop = ifelse(grepl("^A", Sample), "AENP",
                      ifelse(grepl("^B", Sample), "KNP",
                             ifelse(grepl("^H", Sample), "HiP",
                                    ifelse(grepl("^M", Sample), "MNP", "NA")))))

#~~ Define colours
cbPalette <- c("#0066FF", "#CC0000", "#666666", "#FF3300")

#~~ Explicitly state the order of factor variables, so R (ggplot) doesn't use alphabetical order
# For the samples, the order we want to plot them in is, coincidentally, already in alphabetical order,
# so we don't need to explicitly state the order we want, but
# for populations this is not the case (I want KNP to be plotted before HiP). Thus we do the following for
# both psmc_df and boot_df, since we're plotting both of them:
psmc_df$pop <- factor(psmc_df$pop, levels = c("AENP", "KNP", "HiP", "MNP"))
boot_df$pop <- factor(boot_df$pop, levels = c("AENP", "KNP", "HiP", "MNP"))

#~~ Plot per sample, coloured by pop and store as object:
psmc_supp <- ggplot(mapping = aes(dplyr::filter(psmc_df, mu == mu & g == g))) +
  geom_line(aes(YearsAgo, Ne, col = pop), dplyr::filter(psmc_df, mu == mu & g == g), 
            size = 0.7, alpha=0.9) + 
  facet_wrap(~Sample, ncol = 2, scales = "free") +
  geom_line(aes(YearsAgo, Ne, group = Boot, col = pop), boot_df, size = .1, alpha = 0.2) +
  theme_bw() +
  scale_colour_manual(values = cbPalette,
                      name = "Population:",
                      labels = c("AENP", "KNP", "HiP", "MNP")) +
  theme(legend.position = "top") +
  scale_x_log10(label = comma) +
  annotation_logticks(sides = "b", alpha = 0.5) +
  scale_y_continuous(limits = c(-1,2.5e5),
                     label = comma) +
  xlab(expression(paste("Years before present (g = 7.5, ", mu," = 1.5e-8",")"))) + 
  ylab(expression(paste("IICR (scaled in units of 4",italic(N[e]),mu,")")))

# Save as pdf:
pdf(file="figs/psmc_supp.pdf", 
    height=50, width=8)
psmc_supp
dev.off()
