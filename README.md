The project applied comparative genomics and statistical approaches, including phylogenetic reconstruction and pan-genome-wide association analysis, 
to a curated dataset of 107 well-characterised Mycobacterium species to dissect the genetic basis of host adaptation and parasitism.
The bioinformatics pipeline is as follows:\
\
**Dataset and quality control**\
Genomes and metadata were downloaded from ncbi database using the NCBI Datasets command-line tools (CLI) v15.29.0. The dataset was filtered using the quality estimates from CheckM v1.1.2 assessment.
The subsample dataset was manually curated.\
\
**Pan-genome analysis and phylogenetic reconstruction**\
Prokka v1.14.6 was used for the whole genome annotation. Panaroo v1.3.4 was employed in strict mode and a 60% protein family sequence identity threshold for the pan-genome analysis. The ML-tree
was generated and pruned using IQ-TREE v.2.2.6 and ClonalFrameML v1.13. The phylogeny reconstruction was parsed with manually curated metadata and was visualised using R.\
\
**Pangenome-wide association study (panGWAS) analysis**\
Scoary v1.6.16 was employed for the gene-trait association study with 10000 permutations. The significant hits from panGWAS results was analysed and visualised uding R.
