#! /usr/bin/bash
#$ -N checkm
#$ -l h_rt=48:0:0
#$ -l mem=5G
#$ -wd /home/zczl463/Scratch/mb/checkmout
#$ -pe smp 12

source /home/zczl463/anaconda3/etc/profile.d/conda.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/shared/ucl/apps/intel/2022.2/intelpython/python3.9/lib/libimf.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/shared/ucl/apps/miniconda/4.10.3/pkgs/krb5-1.16.4-h2fd8d38_0/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/shared/ucl/apps/miniconda-jhub-test/4.8.3/gnu-4.9.2/lib/
conda activate checkm
export CHECKM_DATA_PATH=/home/zczl463/Scratch/mb/reference
checkm lineage_wf -x fna -t 12 /home/zczl463/Scratch/mb/qc /home/zczl463/Scratch/mb/checkm_output --tab_table -f /home/zczl463/Scratch/mb/checkm_report
