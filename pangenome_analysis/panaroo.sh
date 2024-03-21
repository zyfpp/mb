#!/usr/bin/bash
#$ -N panaroo3
#$ -l h_rt=48:0:0
#$ -l mem=30G
#$ -wd /home/zczl463/Scratch/mb/panout
#$ -pe smp 10

#I have compared a range of protein family thresholds (-f) from 0.30 to 0.95 and used a threshold of 0.60 for further analysis
source /home/zczl463/miniconda3/etc/profile.d/conda.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/shared/ucl/apps/intel/2022.2/intelpython/python3.9/lib/libimf.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/shared/ucl/apps/miniconda/4.10.3/pkgs/krb5-1.16.4-h2fd8d38_0/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/shared/ucl/apps/miniconda-jhub-test/4.8.3/gnu-4.9.2/lib/
conda activate panaroo

panaroo -i /home/zczl463/Scratch/mb/subsetgff/*.gff -o /home/zczl463/Scratch/mb/core_gene/out_3/ --clean-mode strict -a core --aligner mafft --threshold 0.98 -f 0.6 --search_radius 10000 --refind_prop_match 0.8 --core_threshold 0.9 -t 10
