require(dplyr)
require(tidyverse)
require(pheatmap)
require(kableExtra)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(ggrepel)

#read gene-phenotype association files
growth <- read.csv("F:/zyf/mbdata/core_gene_alignment/ML/MLout/scoary_out2/growth_rate_12_03_2024_1911.results.csv")
growth <- growth[, c("Non.unique.Gene.name", "Annotation", "Odds_ratio", "Naive_p", "Benjamini_H_p", "Bonferroni_p", "Empirical_p")]
growth <- growth %>%
  # Remove first semicolons before splitting
  mutate(Non.unique.Gene.name = gsub("^;", "", Non.unique.Gene.name)) %>%
  # Split by the semicolon
  mutate(Non.unique.Gene.name = strsplit(as.character(Non.unique.Gene.name), ";")) %>%
  # create new rows for splitted genes
  unnest(Non.unique.Gene.name) %>%
  # Remove whitespace
  mutate(Non.unique.Gene.name = trimws(Non.unique.Gene.name))
growth <- growth[growth$Non.unique.Gene.name != "", ]


adaptation <- read.csv("F:/zyf/mbdata/core_gene_alignment/ML/MLout/scoary_out2/is_host_adapted_12_03_2024_1911.results.csv")
adaptation <- adaptation[, c("Non.unique.Gene.name", "Annotation", "Odds_ratio", "Naive_p", "Benjamini_H_p", "Bonferroni_p", "Empirical_p")]
adaptation <- adaptation %>%
  mutate(Non.unique.Gene.name = gsub("^;", "", Non.unique.Gene.name)) %>%
  mutate(Non.unique.Gene.name = strsplit(as.character(Non.unique.Gene.name), ";")) %>%
  unnest(Non.unique.Gene.name) %>%
  mutate(Non.unique.Gene.name = trimws(Non.unique.Gene.name))
adaptation <- adaptation[adaptation$Non.unique.Gene.name != "", ]


infection <- read.csv("F:/zyf/mbdata/core_gene_alignment/ML/MLout/scoary_out2/is_human_infection_12_03_2024_1911.results.csv")
infection <- infection[, c("Non.unique.Gene.name", "Annotation", "Odds_ratio", "Naive_p", "Benjamini_H_p", "Bonferroni_p", "Empirical_p")]
infection <- infection %>%
  mutate(Non.unique.Gene.name = gsub("^;", "", Non.unique.Gene.name)) %>%
  mutate(Non.unique.Gene.name = strsplit(as.character(Non.unique.Gene.name), ";")) %>%
  unnest(Non.unique.Gene.name) %>%
  mutate(Non.unique.Gene.name = trimws(Non.unique.Gene.name))
infection <- infection[infection$Non.unique.Gene.name != "", ]


high_temp <- read.csv("F:/zyf/mbdata/core_gene_alignment/ML/MLout/scoary_out2/temp_37_12_03_2024_1911.results.csv")
high_temp <- high_temp[, c("Non.unique.Gene.name", "Annotation", "Odds_ratio", "Naive_p", "Benjamini_H_p", "Bonferroni_p", "Empirical_p")]
high_temp <- high_temp %>%
  mutate(Non.unique.Gene.name = gsub("^;", "", Non.unique.Gene.name)) %>%
  mutate(Non.unique.Gene.name = strsplit(as.character(Non.unique.Gene.name), ";")) %>%
  unnest(Non.unique.Gene.name) %>%
  mutate(Non.unique.Gene.name = trimws(Non.unique.Gene.name))
high_temp <- high_temp[high_temp$Non.unique.Gene.name != "", ]

growth2 <- subset(growth, Empirical_p < 0.05 & Naive_p < 0.05)
adaptation2 <- subset(adaptation, Empirical_p < 0.05 & Naive_p < 0.05)
infection2 <- subset(infection, Empirical_p < 0.05 & Naive_p < 0.05)
high_temp2 <- subset(high_temp, Empirical_p < 0.05 & Naive_p < 0.05)


#Generate a plot to visualise significant association across phenotypes
growth2$Phenotype <- 'Slow growth rate'
infection2$Phenotype <- 'Ability to infect warm-blooded hosts'
adaptation2$Phenotype <- 'Host adaptation'
low_temp2$Phenotype <- 'Grows at 25C'
high_temp2$Phenotype <- 'Grows at 37C'

df <- rbind(growth2, infection2, adaptation2, high_temp2)
df <- df %>%
  mutate(gene_name = gsub("_.*", "", Non.unique.Gene.name))
df <- df %>%
  distinct(gene_name, Phenotype, .keep_all = TRUE)

# Filter to keep only genes appearing in at least three phenotypes
counts <- df %>%
  group_by(gene_name) %>%
  summarise(pheno_count = n_distinct(Phenotype))
abundant_genes <- counts %>%
  filter(pheno_count >= 3) %>%
  select(gene_name)
df <- df %>%
  semi_join(abundant_genes, by = "gene_name")
heatmap_data <- df %>%
  select(gene_name, Phenotype, Odds_ratio) %>%
  pivot_wider(names_from = Phenotype, values_from = Odds_ratio)

gene_names <- heatmap_data$gene_name
heatmap_matrix <- as.matrix(heatmap_data[,-1])
rownames(heatmap_matrix) <- gene_names
heatmap_matrix[is.na(heatmap_matrix)] <- 0
heatmap_matrix <- t(heatmap_matrix)
breaks <- c(0, 0.5, 1, 2, 4, 6, 8, max(heatmap_matrix, na.rm = TRUE) + 0.1)
colours <- c("lightskyblue1", "white", "#FFE3E4", "#FFC1C3", "#FFA3A7", "#FF8288", "#FA4E5A")
heatmap_plot <- pheatmap(heatmap_matrix, clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         clustering_method = "complete",
                         scale = "none", breaks = breaks, color = colours, legend = FALSE,
                         fontsize = 15,       
                         fontsize_row = 18,  
                         fontsize_col = 13)
# Define legends
legend_breaks <- c("0 - 0.5", "0.5 - 1", "1 - 2", "2 - 4", "4 - 6", "6 - 8", ">8")
legend_colors <- c("lightskyblue1", "white", "#FFE3E4", "#FFC1C3", "#FFA3A7", "#FF8288", "#FA4E5A")
plot(0, 0, type = "n", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
legend("topright", legend = legend_breaks, fill = legend_colors, cex = 1.8, title = "Odds ratio")







#identify the 10 most significant hits for each phenotypes
growth3 <- growth2 %>%
  mutate(gene_name = gsub("_.*", "", Non.unique.Gene.name)) %>%
  distinct(gene_name, .keep_all = TRUE)
growth3 <- head(growth3[order(growth3$Empirical_p, decreasing = FALSE), ], 10)

adaptation3 <- adaptation2 %>%
  mutate(gene_name = gsub("_.*", "", Non.unique.Gene.name)) %>%
  distinct(gene_name, .keep_all = TRUE)
adaptation3 <- head(adaptation3[order(adaptation3$Empirical_p, decreasing = FALSE), ], 11)

infection3 <- infection2 %>%
  mutate(gene_name = gsub("_.*", "", Non.unique.Gene.name)) %>%
  distinct(gene_name, .keep_all = TRUE)
infection3 <- head(infection3[order(infection3$Empirical_p, decreasing = FALSE), ], 10)

high_temp3 <- high_temp2 %>%
  mutate(gene_name = gsub("_.*", "", Non.unique.Gene.name)) %>%
  distinct(gene_name, .keep_all = TRUE)
high_temp3 <- head(high_temp3[order(high_temp3$Empirical_p, decreasing = FALSE), ], 10)




#Volcano plot for the panGWAS analysis
#annotate top 5 hits for each phenotype with odds ratio greater or smaller than 1
growth2$source <- 'growth'
infection2$source <- 'infection'
adaptation2$source <- 'adaptation'
high_temp2$source <- 'high_temp'
combined_df <- rbind(growth2, infection2, adaptation2, high_temp2)
combined_df$neg_log10_empirical_p <- -log10(combined_df$Empirical_p) #transform p to negative log value
combined_df$log_Odds_ratio <- log(combined_df$Odds_ratio) #transform odds ratio to log value
combined_df <- combined_df[combined_df$Odds_ratio != 'Inf' & combined_df$Odds_ratio > 0, ] #remove extremely large and small value for visualisation

#create labels for top 5 hits
acquisition <- combined_df %>%
  filter(Odds_ratio > 1) %>%
  arrange(Phenotype, neg_log10_empirical_p) %>%
  group_by(Phenotype) %>%
  slice_max(neg_log10_empirical_p, n = 5) 
deletion <- combined_df %>%
  filter(Odds_ratio < 1) %>%
  arrange(Phenotype, neg_log10_empirical_p) %>%
  group_by(Phenotype) %>%
  slice_max(neg_log10_empirical_p, n = 5)
annotations <- rbind(acquisition, deletion)

p <- ggplot(data=combined_df, aes(x=log_Odds_ratio, y=neg_log10_empirical_p, color=Phenotype)) +
  geom_point() +
  theme_minimal() +
  geom_vline(xintercept = log(1), col = "blue", linetype="dashed") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype="dashed") +
  geom_label_repel(data = annotations, 
                   aes(x=log_Odds_ratio, y=neg_log10_empirical_p, label = gene_name),
                   size = 8,  # Ensure this is defined or replace with a numeric value
                   box.padding = 0.35, 
                   point.padding = 0.4, 
                   min.segment.length = 0.1,
                   segment.color = 'grey50',
                   force = 10, 
                   max.overlaps = Inf,
                   show.legend = FALSE,
                   fill = 'white') +  # Background box color
  scale_color_manual(values = c('growth' = 'seagreen2', 'infection' = 'firebrick', 
                                'adaptation' = 'dodgerblue4', 'high_temp' = 'mediumpurple1')) +
  labs(x = "Log(Odds Ratio)", y = "-log10(Empirical P-Value)")

p <- p + theme(legend.title = element_text(size = 26),
               axis.text = element_text(size = 20),    
               axis.title = element_text(size = 20),
               legend.text = element_text(size = 24))
print(p)









