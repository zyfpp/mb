#the ML-tree was generated from core gene alignments using IQtree 
#the rooted tree was pruned to account for recombination using ClonalFrameML
iqtree2 -s core_gene_alignment.masked_to_98525pos.aln -m MFP -nt 10 --prefix f60
ClonalFrameML rooted.tree core_gene_alignment.masked_to_98525pos.aln recombination_pruned 
