#BiocManager::install("ggtreeExtra")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(ggtreeExtra)
require(ggnewscale)
require(ggutils)
require(randomcoloR)
require(foreach)
require(Biostrings)
require(doParallel)

#alignment QC
rm(list = ls())
setwd("F:/zyf/mbdata/core_gene_alignment/ML/MLout")
aln_path <- "core_gene_alignment.aln"
aln <- readDNAStringSet(aln_path)
aln_mat <- as.matrix(aln)

# Get gapped positions
gappy_pos <- apply(aln_mat, 2,
                   function(x) {sum(x %in% c("-", "N")) / length(x) > 0.1})
n_unmasked <- sum(!gappy_pos)

# Mask gapped positions
aln_mat[, gappy_pos] <- "N"
ncol(aln_mat) == unique(width(aln))
masked_aln <- DNAStringSet(apply(aln_mat, 1, paste0, collapse = ""))

# Save masked alignment
save_name <- gsub("\\.aln", str_glue("\\.masked_to_{n_unmasked}pos.aln"), aln_path)
writeXStringSet(masked_aln, save_name)

meta <- read.csv("sub_1genome.csv", check.names = F, stringsAsFactors = F)
tree <- read.tree("recombination_pruned.labelled_tree.newick")



#the ML-tree was generated from core gene alignments using IQtree 
#the rooted tree was pruned to account for recombination using ClonalFrameML
#root to M.grossiae
outgrp <- meta.match %>% 
  filter(spc == "grossiae")
outgrp <- c(outgrp$genome)
tree$tip.label <- sub("^(.*?_.*?)_.*$", "\\1", tree$tip.label)
rooted <- root(tree, outgroup = outgrp, resolve.root = T)
write.tree(rooted, "F:/zyf/mbdata/core_gene_alignment/rooted.tree")

#Match metadata to tips
meta.match <- tibble(genome = tree$tip.label) %>%
  left_join(meta.match) %>%
  mutate(tip_label = genome) %>%
  mutate(host = ifelse(host == "", NA, host),
         spc = ifelse(spc == "", NA, spc))

all(tree$tip.label == meta.match$genome)

#Plot tree
col_pal <- c("seagreen4", "darkgoldenrod3","#d1495b", "lightsalmon", 
             "maroon4", "darkorchid4", "slateblue4",
             "lightpink4", "khaki", "lightpink", 
             "navajowhite3", "darkolivegreen4", "lightgreen",
             "tan2", "blue", "firebrick4")
col_pal2 <- c("darkgrey","#2e4057")
col_pal3 <- c("cornsilk3", "darkcyan")
col_pal4 <- c("lightsteelblue", "mediumpurple")
col_pal5 <- c("dodgerblue1","dodgerblue4")
set.seed(66)

p <- ggtree(tree,
            layout = 'circular', 
            size= 0.1,
            color = "darkslategrey",
            options(ignore.negative.edge = TRUE)) %<+% meta.match +
  geom_tippoint(aes(hjust = 0.5, color = growth_rate), 
                alpha = 1, 
                size = 1)

p <- p + scale_color_discrete(na.translate = F, guide = guide_legend(order = 1))

#host layer
p <- p + geom_fruit(geom = geom_tile,
                    aes(fill = host),
                    offset = 0.6,
                    width = 0.01) +
  scale_fill_manual(values = col_pal, 
                    name = "Host",
                    na.translate = F,
                    guide = guide_legend(position = "bottom")) #order = 2

#adaptation layer
p <- p + new_scale_fill() +
  geom_fruit(geom = geom_tile,
             aes(fill = is_host_adapted),
             offset = 0.15,
             width = 0.01) +
  scale_fill_manual(values = col_pal2, 
                    name = "Adaptation",
                    na.translate = F,
                    guide = guide_legend(order = 3))
#infection layer
p <- p + new_scale_fill() +
  geom_fruit(geom = geom_tile,
             aes(fill = is_human_infection),
             offset = 0.15,
             width = 0.01) +
  scale_fill_manual(values = col_pal3, 
                    name = "Infects warm-blooded hosts?",
                    guide = guide_legend(order = 4))

#temp layers
p <- p + new_scale_fill() +
  geom_fruit(geom = geom_tile,
             aes(fill = temp_25),
             offset = 0.15,
             width = 0.01) +
  scale_fill_manual(values = col_pal4, 
                    name = "Grows at 25C?",
                    na.translate = F,
                    guide = guide_legend(order = 5))

p <- p + new_scale_fill() +
  geom_fruit(geom = geom_tile,
             aes(fill = temp_37),
             offset = 0.15,
             width = 0.01) +
  scale_fill_manual(values = col_pal5, 
                    name = "Grows at 37C?",
                    na.translate = F,
                    guide = guide_legend(order = 6))

# Add labels
p <- p + labs(color = "Growth rate")
p <- p + geom_tiplab(aes(label = spc), align = TRUE, linetype = "solid", size = 2)
p <- p + theme(legend.title = element_text(size = 20),
               legend.text = element_text(size = 16))
print(p)
