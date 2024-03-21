library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork)
library(ggrepel)

#plot species-phenotype distribution of my subsample
meta.match <- read.csv("F:/zyf/mbdata/core_gene_alignment/ML/MLout/matched_meta.csv")
text_size = 6 
title_text_size = 22 
legend_text_size = 18 

#Distribution of species by growth rate
p1 <- ggplot(df1, aes(x = "", y = distinct_spc, fill = growth_rate)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  theme_void() +
  scale_fill_brewer(palette = "Set1") +
  labs(fill = "Growth Rate", title = "Distribution of species by growth rate") +
  geom_text_repel(aes(x = 1, y = cumsum(distinct_spc) - distinct_spc / 2,
                      label = paste0(distinct_spc, " (", round(distinct_spc / sum(distinct_spc) * 100, 1), "%)")),
                  nudge_x = 0.15, size = text_size) +
  theme(legend.title = element_text(size = legend_text_size), 
        legend.text = element_text(size = legend_text_size),
        plot.title = element_text(size = title_text_size))

#Distribution of species by host adaptation
p2 <- ggplot(df2, aes(x = "", y = distinct_spc, fill = is_host_adapted)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  theme_void() +
  scale_fill_brewer(palette = "Set2") +
  labs(fill = "Host adaptation", title = "Distribution of species by host adaptation") +
  geom_text_repel(aes(label = paste0(distinct_spc, " (", round(distinct_spc / sum(distinct_spc) * 100, 1), "%)")),
                  position = position_stack(vjust = 0.5), size = text_size) +
  theme(legend.title = element_text(size = legend_text_size), 
        legend.text = element_text(size = legend_text_size),
        plot.title = element_text(size = title_text_size))

#Distribution of species by infection
p3 <- ggplot(df3, aes(x = "", y = distinct_spc, fill = is_human_infection)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  theme_void() +
  scale_fill_brewer(palette = "Set3") +
  labs(fill = "Reported infection in warm-blooded host?", title = "Distribution of species by infection") +
  geom_text_repel(aes(label = paste0(distinct_spc, " (", round(distinct_spc / sum(distinct_spc) * 100, 1), "%)")),
                  position = position_stack(vjust = 0.5), size = text_size) +
  theme(legend.title = element_text(size = legend_text_size), 
        legend.text = element_text(size = legend_text_size),
        plot.title = element_text(size = title_text_size))

#Distribution of species by growth at 37C
p4 <- ggplot(df4, aes(x = "", y = distinct_spc, fill = temp_37)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  theme_void() +
  scale_fill_brewer(palette = "Set4") +
  labs(fill = "Grows at 37C?", title = "Distribution of species by growth at 37C") +
  geom_text_repel(aes(label = paste0(distinct_spc, " (", round(distinct_spc / sum(distinct_spc) * 100, 1), "%)")),
                  position = position_stack(vjust = 0.5), size = text_size) +
  theme(legend.title = element_text(size = legend_text_size), 
        legend.text = element_text(size = legend_text_size),
        plot.title = element_text(size = title_text_size))

# Combine and print the plots
combined_plot <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
print(combined_plot)