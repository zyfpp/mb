#pan-GWAS with 10000 permutations
#the gene-presence/absence file was from the pan-genome reesults
#a list of phenotypes was manually coded binarily from metadata

scoary -g gene_presence_absence_roary.csv -t traits.csv -o scoary_out -e 10000 --threads 12
