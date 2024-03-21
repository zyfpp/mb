#BiocManager::install("ggtreeExtra")
require(tidyverse)
require(ggplot2)
require(ggtreeExtra)
require(ggnewscale)
require(ggutils)



#assessing panaroo protein family threshold parameters
#the plot identified a threshold of 0.60 to discriminate between core (>95% prevalence across species) and accessory (<95%) genes. 
scaling_factor <- mean(core_accessory$accessory_gene) / mean(core_accessory$core_gene) / 5
core_accessory$core_gene_scaled <- core_accessory$core_gene * scaling_factor
df <- melt(core_accessory, id.vars = "P_family_threshold",
                            measure.vars = c("core_gene_scaled", "accessory_gene"))

df$variable <- sub("core_gene_scaled", "core_gene", df$variable)


p <- ggplot(df, aes(x = as.factor(P_family_threshold), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  geom_text(aes(label = ifelse(variable == "core_gene", round(core_accessory$core_gene), round(value))), 
            position = position_dodge(width = 0.7), vjust = -0.5, size = 6) +
  scale_fill_manual(values = c("core_gene" = "cornflowerblue", "accessory_gene" = "cornsilk3")) +
  labs(x = "Protein Family Threshold", y = "Gene Clusters Count") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 14),  # Increase x-axis numbers size
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  theme(legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20),
        plot.title = element_text(size = 20))

print(p)