# This script filters the given dataset to produce a summarized visualization
# of key variables within the dataset.
#
# Zhuangzhuang Geng, D3B 2024

## load libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggthemes)

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Define the file.path to output directories
analysis_dir <- file.path(root_dir, "figures", "manuscript_OPC")
plots_dir <- file.path(analysis_dir, "plots")

## source theme_Publication
source(file.path(analysis_dir, "utils", "theme_for_plot.R"))

# Create directories to hold the output.
if (!dir.exists(analysis_dir)) {
  dir.create(analysis_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Read in dataset and remove NAs
histologies_df <- readr::read_tsv(file.path(root_dir, "data", "histologies.tsv"), 
                                  guess_max = 10000) %>%
  # mutate(group = case_when(cohort %in% c("TCGA", "GTEx") ~ "Adult", TRUE ~ "Pediatric")) %>%
  mutate(group = case_when(cohort == "TCGA" ~ "Adult tumors",
                           cohort == "GTEx" ~ "Adult normal", 
                           sub_cohort == "CPTAC GBM" ~ "Adult tumors",
                           !is.na(pathology_diagnosis) & !cohort %in% c("TCGA", "GTEx") ~ "Pediatric",
                           TRUE ~ NA_character_))

## Generate the pie plot for tumor description
tumor_descriptor_ped <- histologies_df %>%
  filter(group == "Pediatric") %>%
  mutate(tumor_descriptor = case_when(tumor_descriptor == "Initial CNS Tumor" ~ "Primary Tumor",
                                      tumor_descriptor %in% c("Recurrence", "Relapse") ~ "Recurrence or Relapse",
                                      tumor_descriptor == "Progressive Disease Post-Mortem" ~ "Progressive",
                                      TRUE ~ tumor_descriptor)) %>%
  select(match_id, tumor_descriptor) %>% 
  distinct() %>% 
  group_by(tumor_descriptor) %>% 
  tally() %>% 
  filter(!(tumor_descriptor %in% c("Unknown", NA))) %>% 
  mutate(percentage = paste0(round((n/sum(n))*100, digit = 2), "%")) %>%
  arrange(desc(n))

tumor_descriptor_adult <- histologies_df %>%
  filter(cohort == "TCGA") %>%
  mutate(tumor_descriptor = case_when(tumor_descriptor == "Initial CNS Tumor" ~ "Primary Tumor",
                                      tumor_descriptor %in% c("Recurrence", "Relapse") ~ "Recurrence or Relapse",
                                      tumor_descriptor == "Progressive Disease Post-Mortem" ~ "Progressive",
                                      TRUE ~ tumor_descriptor)) %>%
  group_by(tumor_descriptor) %>% 
  tally() %>% 
  filter(!(tumor_descriptor %in% c("Unknown", NA))) %>% 
  mutate(percentage = paste0(round((n/sum(n))*100, digit = 2), "%")) %>%
  arrange(desc(n))

tum_colors <- c("Primary Tumor" = "blueviolet", 
                "Recurrence or Relapse" = "darkslategray2",
                "Progressive" = "darkseagreen",
                "Deceased" = "deeppink",
                "Second malignancy" = "salmon", 
                "Metastatic" = "darkgoldenrod2", 
                "Post-treatment" = "yellow1", 
                "Residual" = "magenta2")

ped_tum_plot <- ggplot(tumor_descriptor_ped, aes(x = "", y = n, fill = factor(tumor_descriptor, levels = tumor_descriptor))) + 
  geom_bar(stat = "identity", width = 1, color = "white") + 
  coord_polar("y", start = 0) + 
  scale_fill_manual(name = "Tumor Event", 
                    labels = paste0(tumor_descriptor_ped$tumor_descriptor, " (", tumor_descriptor_ped$percentage, ")"), 
                    values = tum_colors) +
  theme_void() + 
  theme(legend.title = element_text(face = "bold", size = 17),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold")) + # Customize legend title text if needed
  labs(title = "Pediatric Tumors")

tcga_plot <- ggplot(tumor_descriptor_adult, aes(x = "", y = n, fill = factor(tumor_descriptor, levels = tumor_descriptor))) + 
  geom_bar(stat = "identity", width = 1, color = "white") + 
  coord_polar("y", start = 0) + 
  scale_fill_manual(name = "Tumor Event", 
                    labels = paste0(tumor_descriptor_adult$tumor_descriptor, " (", tumor_descriptor_adult$percentage, ")"), 
                    values = tum_colors) + 
  theme_void() + 
  theme(legend.title = element_text(face = "bold", size = 17),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold")) +  # Customize legend title text if needed
  labs(title = "Adult tumors")


# Combine the plots into one
combined_td_plot <- plot_grid(ped_tum_plot, tcga_plot, 
                              ncol = 1, align = "v", rel_heights = c(2, 2)  # Adjust relative heights
)

# Display the combined plot
print(combined_td_plot)

ggsave(file.path(plots_dir, "tumor_descriptor.png"), height = 8, width = 10)

## generate stacked barplot for experimental strategy
exp_strategy_df <- histologies_df %>%
  filter(!is.na(group)) %>%
  #mutate(group = factor(group, levels = c("Pediatric", "Adult"))) %>%
  group_by(match_id) %>% 
  select(group, cohort, experimental_strategy)

exp_level <- unique(exp_strategy_df$experimental_strategy)[order(unique(exp_strategy_df$experimental_strategy))]

# table to be plotted
table(exp_strategy_df$cohort, exp_strategy_df$group)

# Create the plot for the pediatric cohort
plot_pediatric <- ggplot(subset(exp_strategy_df, group == "Pediatric"), aes(x = cohort, fill = factor(experimental_strategy, levels = exp_level))) + 
  geom_bar() +  # Default width
  scale_fill_manual("Experimental Strategy", values = c("Methylation" = "darkgoldenrod2", 
                                                        "miRNA-Seq" = "yellow1", 
                                                        "Phospho-Proteomics" = "magenta2",
                                                        "RNA-Seq" = "salmon", 
                                                        "Targeted Sequencing" = "darkslategray2",
                                                        "WGS" = "blueviolet", 
                                                       "Whole Cell Proteomics"= "deeppink",
                                                       "WXS" = "darkseagreen")) +
  coord_flip() +
  ylim(0,10000) +
  theme_Publication() +
  labs(title = "Pediatric Tumors", x = "", y = "Number of tumors")  # Add axis labels


# Create the plot for the adult tumors
plot_adult <- ggplot(subset(exp_strategy_df, group == "Adult tumors"), aes(x = cohort, fill = factor(experimental_strategy, levels = exp_level))) + 
  geom_bar(width = 0.5) +  # Adjust the width here
  scale_fill_manual("Experimental Strategy", values = c("Methylation" = "darkgoldenrod2", 
                                                        "miRNA-Seq" = "yellow1", 
                                                        "Phospho-Proteomics" = "magenta2",
                                                        "RNA-Seq" = "salmon", 
                                                        "Targeted Sequencing" = "darkslategray2",
                                                        "WGS" = "blueviolet", 
                                                        "Whole Cell Proteomics"= "deeppink",
                                                        "WXS" = "darkseagreen")) +
  coord_flip() +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
  labs(title = "Adult Tumors", x = "", y = "Number of samples")  # Add axis labels


# Create the plot for the adult cohort
plot_adult2 <- ggplot(subset(exp_strategy_df, group == "Adult normal"), aes(x = cohort, fill = factor(experimental_strategy, levels = exp_level))) + 
  geom_bar(width = 0.5) +  # Adjust the width here
  scale_fill_manual("Experimental Strategy", values = c("Methylation" = "darkgoldenrod2", 
                                                        "miRNA-Seq" = "yellow1", 
                                                        "Phospho-Proteomics" = "magenta2",
                                                        "RNA-Seq" = "salmon", 
                                                        "Targeted Sequencing" = "darkslategray2",
                                                        "WGS" = "blueviolet", 
                                                        "Whole Cell Proteomics"= "deeppink",
                                                        "WXS" = "darkseagreen")) +
  coord_flip() +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
  labs(title = "Adult normal tissues", x = "", y = "Number of samples")  # Add axis labels

# Combine the plots into one
combined_plot <- plot_grid(plot_pediatric, plot_adult, plot_adult2, 
                           ncol = 1, align = "v", rel_heights = c(2, 1.3, 1.1)  # Adjust relative heights
)

# Display the combined plot
print(combined_plot)


ggsave(file.path(plots_dir, "experimental_strategy.png"), height = 8, width = 6.5)

